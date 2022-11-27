# -*- coding: utf-8 -*-
"""
Written by: Yasuhiro Miyazawa
"""

# Import packages
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt

from scipy import optimize
from scipy import interpolate

import sys
import os
import shutil
import zipfile
import glob

from tqdm import tqdm

# Figure parameters
plt.style.use('./common/custom.mplstyle')


class StartuckAnalysisTool:
    def __init__(self,):
        print('{0:*^49}'.format(''))
        print('* {0:^45} *'.format('Startuck Analysis Tool'))
        print('{0:*^49}'.format(''))

        pass

    def calc_rotationMatrix_X(self, th):
        rotT = np.array([[1, 0, 0], [0, np.cos(th), -np.sin(th), 0], [0, np.sin(th), np.cos(th)]])

        return rotT

    def calc_rotationMatrix_Y(self, th):
        rotT = np.array([[np.cos(th), 0, -np.sin(th), 0], [0, 1, 0], [np.sin(th), 0, np.cos(th)]])

        return rotT

    def calc_rotationMatrix_Z(self, th):
        rotT = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])

        return rotT

    def calc_rotationMatrix_rodrig(self, vec, th):
        vec_k = vec / np.linalg.norm(vec)
        mat_K = np.array([[0, -vec_k[2], vec_k[1]], [vec_k[2], 0, -vec_k[0]], [-vec_k[1], vec_k[0], 0]])
        rotRodrig = np.identity(3) + np.sin(th) * mat_K + (1. - np.cos(th)) * (mat_K @ mat_K)

        return rotRodrig

    def geo_init_startuck(self, La: float,
                          theta_M0: float,
                          st_type: str = 'tri',
                          fig_out: bool = True,
                          ):
        # **********************************************
        #        Startuck unit node definition
        # **********************************************

        self.La = La
        self.theta_M0 = theta_M0
        match st_type:
            case 'tri':
                st_profile = 'Triangle'
                n_nodepoly = 3
                self.n_node = 7
                self.n_edge = 12
                self.n_poly = 6

                self.Lb = np.sqrt(3.) * La / 3.
                self.Lc = 2 * np.sqrt(3.) * La / 3.

                # Define edge connection
                EdgeConct = np.zeros((self.n_edge, 2), dtype=int)
                EdgeConct[0, :] = [0, 1]
                EdgeConct[1, :] = [0, 2]
                EdgeConct[2, :] = [0, 3]
                EdgeConct[3, :] = [0, 4]
                EdgeConct[4, :] = [0, 5]
                EdgeConct[5, :] = [0, 6]
                EdgeConct[6, :] = [1, 2]
                EdgeConct[7, :] = [2, 3]
                EdgeConct[8, :] = [3, 4]
                EdgeConct[9, :] = [4, 5]
                EdgeConct[10, :] = [5, 6]
                EdgeConct[11, :] = [6, 1]

                # Define polygon
                Polyg = np.zeros((self.n_poly, n_nodepoly), dtype=int)
                Polyg[0, :] = [0, 1, 2]
                Polyg[1, :] = [0, 2, 3]
                Polyg[2, :] = [0, 3, 4]
                Polyg[3, :] = [0, 4, 5]
                Polyg[4, :] = [0, 5, 6]
                Polyg[5, :] = [0, 6, 1]

                # Useful rotation matrices
                self.rotTz60 = self.calc_rotationMatrix_Z(np.pi / 3)
                self.rotTz120 = self.calc_rotationMatrix_Z(2 * np.pi / 3)
                self.rotTz60n = self.calc_rotationMatrix_Z(-np.pi / 3)
                self.rotTz120n = self.calc_rotationMatrix_Z(-2 * np.pi / 3)

                self.find_geo = self.find_geo_triangle

        # Define nodal coordinates
        vert_xyz = self.find_geo(theta_M0)

        self.vert_xyz0 = vert_xyz
        self.EdgeConct = EdgeConct
        self.Polyg = Polyg
        self.st_profile = st_profile

        print('{0:*^49}'.format(''))
        print(' Initialize Startuck Fold Geometry')
        print('  {0:<24s} : {1:24s}'.format('Startuck profile', st_profile))
        print('  {0:<24s} : {1:.6f}'.format('Side length', self.La))
        print('  {0:<24s} : {1:d} {2:d} {3:d}'.format('Node, crease, polygon', self.n_node, self.n_edge, self.n_poly))
        print('  {0:<24s} : {1:.6f} (deg)'.format('Initial fold angle', np.degrees(self.theta_M0)))
        print('{0:*^49}'.format(''))

        if fig_out:
            from .plot_geometry import plot_projection
            plot_projection(vert_xyz=vert_xyz,
                            Polyg=self.Polyg,
                            EdgeConct=self.EdgeConct,
                            figname='(Startuck-init)')

        return

    def calc_thetaS(self, theta_M: float | npt.ArrayLike):
        thGS = np.arccos(np.sin(theta_M))
        thHS = 0.5 * np.arccos(np.cos(2 * theta_M))
        thHH = np.emath.arcsin(np.sin(theta_M) / (np.sqrt(4. - 3. * np.sin(theta_M)**2))).real
        thHM = np.emath.arccos(1. / (np.sqrt(4. - 3. * np.sin(theta_M)**2))).real
        psi = np.pi - thGS - thHS - thHH - thHM
        theta_S = np.arcsin(np.cos(psi))

        return theta_S

    def calc_thetaG(self, theta_M: float | npt.ArrayLike):
        theta_S = self.calc_thetaS(theta_M)
        theta_GM = np.emath.arccos(np.sin(theta_S))
        theta_GS = np.emath.arccos(np.sin(theta_M))

        return theta_GM, theta_GS

    def calc_thetaH(self, theta_M: float | npt.ArrayLike):
        theta_HM = np.emath.arccos(1. / (np.sqrt(4. - 3. * np.sin(theta_M)**2))).real
        theta_HS = 0.5 * np.arccos(np.cos(2 * theta_M))
        theta_HH = np.emath.arcsin(np.sin(theta_M) / (np.sqrt(4. - 3. * np.sin(theta_M)**2))).real

        return theta_HM, theta_HS, theta_HH

    def calc_evectors(self, theta_M):
        theta_GM, theta_GS = self.calc_thetaG(theta_M)
        e2 = np.array([0.0, np.cos(theta_GM), np.sin(theta_GM)])
        e5 = np.array([0.0, -np.cos(theta_GS), np.sin(theta_GS)])
        e1 = self.rotTz120 @ e5
        e3 = self.rotTz120 @ e1
        e4 = self.rotTz120 @ e2
        e6 = self.rotTz120 @ e4

        return e1, e2, e3, e4, e5, e6

    def find_geo_triangle(self, theta_M: float):

        vert_xyz = np.zeros((self.n_node, 3))
        e1, e2, e3, e4, e5, e6 = self.calc_evectors(theta_M)
        vert_xyz[0, :] = 0.0
        vert_xyz[1, :] = self.Lb * e1
        vert_xyz[2, :] = self.Lc * e2
        vert_xyz[3, :] = self.Lb * e3
        vert_xyz[4, :] = self.Lc * e4
        vert_xyz[5, :] = self.Lb * e5
        vert_xyz[6, :] = self.Lc * e6

        return vert_xyz

    def calc_facet_area(self, vert_xyz: npt.ArrayLike):
        facet_area = np.zeros(self.n_poly)
        for ip in range(self.n_poly):
            vec1 = vert_xyz[self.Polyg[ip, 0], :] - vert_xyz[self.Polyg[ip, 1], :]
            vec2 = vert_xyz[self.Polyg[ip, 2], :] - vert_xyz[self.Polyg[ip, 1], :]
            facet_area[ip] = np.linalg.norm(np.cross(vec1, vec2))

        return facet_area

    def write_vtk(self, fnum: int, vert_xyz: npt.ArrayLike, strain: npt.ArrayLike, st_type: str = 'mono'):
        '''
        write_vtk
            Exports vtk file of the structure.
            Created vtk files can be imported into, e.g., Paraview (https://www.paraview.org/download/).

        Args:
            fnum (int): File number in integer.
            vert_xyz (npt.ArrayLike): xyz coordinates of all nodes
            strain (npt.ArrayLike): Value used for coloring the polygons
        '''

        # Open file
        with open('%s/%s_%05d.vtk' % (self.dir_save, self.fname_vtk, fnum), 'w') as f:
            # Write header lines
            f.write('# vtk DataFile Version 3.0\n')
            f.write('%s Startuck\n' % self.st_profile)
            f.write('ASCII\n')
            num_points = np.size(vert_xyz, axis=0)
            f.write('DATASET POLYDATA\n')
            # Nodes
            # Specify number of nodes
            f.write('POINTS %d double\n' % num_points)
            # Write x,y,z coodinates of all nodes
            for iv in range(num_points):
                f.write("%.15e %.15e %.15e\n" % (vert_xyz[iv, 0], vert_xyz[iv, 1], vert_xyz[iv, 2]))

            #   POLYGONS
            num_dataset = self.n_poly
            num_datanum = 4 * num_dataset
            f.write('POLYGONS %d %d\n' % (num_dataset, num_datanum))
            for ip in range(self.n_poly):
                f.write('3 %d %d %d\n' % tuple(self.Polyg[ip, :]))

            f.write('CELL_DATA %d\n' % num_dataset)
            f.write('SCALARS cell_scalars double\n')
            f.write('LOOKUP_TABLE default\n')
            for i in range(self.n_poly):
                f.write('%.15e\n' % strain[i])

        return

    def create_3Dmodel(self, theta_M: list | np.ndarray,
                       save_zip: bool = False):

        niter = len(theta_M)
        # VTK export
        self.dir_save = './vtk_startuck_{:s}'.format(self.st_profile)
        self.fname_vtk = 'startuck_{:s}'.format(self.st_profile)

        print('Create VTK...')
        # Delete previous data
        if os.path.exists(self.dir_save):
            shutil.rmtree(self.dir_save)
        os.makedirs(self.dir_save)
        # Generate vtk file
        for ii in tqdm(range(niter)):
            vert_xyz = self.find_geo_triangle(theta_M[ii])
            # Dummy array for panel color
            strain = self.calc_facet_area(vert_xyz=vert_xyz)
            self.write_vtk(ii, vert_xyz, strain)

        if save_zip:
            zp = zipfile.ZipFile('%s.zip' % self.dir_save, 'w')
            dfile = glob.glob('%s/*.vtk' % self.dir_save)
            dfile = np.sort(dfile)
            for i in range(len(dfile)):
                zp.write(filename=dfile[i], arcname=None, compress_type=None, compresslevel=9)
            zp.close()

        return


if __name__ == "__main__":

    sat = StartuckAnalysisTool()
    nbin = 256
    La = 1.0
    theta_M0 = np.pi / 2.
    theta_M1 = 0.0
    theta_M = np.linspace(theta_M0, theta_M1, nbin)
    sat.geo_init_startuck(La, theta_M0, st_type='tri')

    theta_S = sat.calc_thetaS(theta_M)

    # plt.figure('theta_MS', figsize=(8, 3))
    # plt.plot(np.degrees(theta_M), np.degrees(theta_S))
    # plt.xlabel('$\\theta_M$ (deg)')
    # plt.ylabel('$\\theta_S$ (deg)')
    # plt.tight_layout()
    # plt.show()
    # exit()

    sat.create_3Dmodel(theta_M=theta_M, save_zip=False)
    sat.create_3Dmodel_unit(theta_M=theta_M, save_zip=False)
