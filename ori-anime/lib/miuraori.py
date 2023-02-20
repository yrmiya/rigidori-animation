# -*- coding: utf-8 -*-
"""
Written by: Yasuhiro Miyazawa
"""

# Import packages
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt

import sys
import os
import shutil
import zipfile
import glob

from tqdm import tqdm

# Figure parameters
# plt.style.use('./common/custom.mplstyle')


class MiuraOriAnalysis:
    def __init__(self, dir_save, fname_vtk):
        print('{0:*^49}'.format(''))
        print('* {0:^45} *'.format('Miura-ori Static Analysis Tool'))
        print('{0:*^49}'.format(''))

        # VTK export
        self.dir_save = dir_save
        self.fname_vtk = fname_vtk

        pass

    def geo_init_miura(self, n_unitx: int, n_unity: int,
                       La: float, Lb: float, Ld: float,
                       alpha: float,
                       theta_M0: float,
                       fig_out: bool = True,
                       ):
        # **********************************************
        #        Miura-ori unit node definition
        # **********************************************
        #  O ----> y
        #  |
        #  | x
        #  v
        #          0 ------------- 3 ---------- 6
        #         /               /            /
        #        /               /            /
        #       /               /            /
        #      1 -------a----- 4 -----b---- 7      ---
        #       \               \            \      |
        #        \               \            \     | d
        #         \               \      alpha \    |
        #          2 ------------- 5 ---------- 8  ---
        # **********************************************

        self.La = La
        self.Lb = Lb
        self.Ld = Ld
        self.alpha = alpha
        self.theta_M0 = theta_M0
        self.n_node = (2 * n_unitx + 3) * n_unity
        self.n_edge = (2 * n_unitx + 1) * (n_unity + 1) + (2 * n_unity + 1) * (n_unitx + 1)
        self.n_poly = 4 * n_unitx * n_unity

        # Define nodal coordinates
        vert_xyz = self.find_geo_miura_unit(theta_M0)

        # Define edge connection
        EdgeConct = np.zeros((self.n_edge, 2), dtype=int)
        EdgeConct[0, :] = [0, 1]
        EdgeConct[1, :] = [1, 2]
        EdgeConct[2, :] = [3, 4]
        EdgeConct[3, :] = [4, 5]
        EdgeConct[4, :] = [6, 7]
        EdgeConct[5, :] = [7, 8]
        EdgeConct[6, :] = [0, 3]
        EdgeConct[7, :] = [3, 6]
        EdgeConct[8, :] = [1, 4]
        EdgeConct[9, :] = [4, 7]
        EdgeConct[10, :] = [2, 5]
        EdgeConct[11, :] = [5, 8]

        # Define polygon
        Polyg = np.zeros((self.n_poly, 4), dtype=int)
        Polyg[0, :] = [0, 1, 4, 3]
        Polyg[1, :] = [1, 2, 5, 4]
        Polyg[2, :] = [3, 4, 7, 6]
        Polyg[3, :] = [4, 5, 8, 7]

        self.vert_xyz0 = vert_xyz
        self.EdgeConct = EdgeConct
        self.Polyg = Polyg

        print('{0:*^49}'.format(''))
        print(' Initialize Miura-ori Fold Geometry')
        print('  {0:<24s} : {1:.6f}'.format('Side length', self.La))
        print('  {0:<24s} : {1:d} {2:d} {3:d}'.format('Node, crease, polygon', self.n_node, self.n_edge, self.n_poly))
        print('  {0:<24s} : {1:.6f} (deg)'.format('Initial fold angle', np.degrees(self.theta_M0)))
        print('{0:*^49}'.format(''))

        if fig_out:
            from .plot_geometry import plot_projection
            plot_projection(vert_xyz=vert_xyz,
                            Polyg=self.Polyg,
                            EdgeConct=self.EdgeConct,
                            figname='(Miura-init)')

        return

    def find_geo_miura_unit(self, theta_M: float):

        vert_xyz = np.zeros((9, 3))
        xi = np.tan(self.alpha) * np.cos(theta_M)
        cota = 1. / np.tan(self.alpha)
        ex = np.array([1., 0, 0.])
        ey = np.array([0., 1, 0.])
        ez = np.array([0., 0, 1.])
        e1 = np.array([0, (1 - xi**2) / (1 + xi**2), (2 * xi) / (1 + xi**2)])
        e2 = 1. / (1. + cota**2) * np.array([np.sin(theta_M), cota, np.cos(theta_M)])
        e2n = 1. / (1. + cota**2) * np.array([-np.sin(theta_M), cota, np.cos(theta_M)])
        vert_xyz[0, :] = self.Ld / np.sin(self.alpha) * e2n
        vert_xyz[1, :] = 0.0
        vert_xyz[2, :] = self.Ld / np.sin(self.alpha) * e2
        vert_xyz[3, :] = self.La * ey + self.Ld / np.sin(self.alpha) * e2n
        vert_xyz[4, :] = self.La * ey
        vert_xyz[5, :] = self.La * ey + self.Ld / np.sin(self.alpha) * e2
        vert_xyz[6, :] = self.La * ey + self.Lb * e1 + self.Ld / np.sin(self.alpha) * e2n
        vert_xyz[7, :] = self.La * ey + self.Lb * e1
        vert_xyz[8, :] = self.La * ey + self.Lb * e1 + self.Ld / np.sin(self.alpha) * e2

        return vert_xyz

    def write_vtk(self, fnum: int, vert_xyz: npt.ArrayLike, strain: npt.ArrayLike):
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
            f.write('Miura ori\n')
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
            num_datanum = 5 * num_dataset
            f.write('POLYGONS %d %d\n' % (num_dataset, num_datanum))
            for ip in range(self.n_poly):
                f.write('4 %d %d %d %d\n' % tuple(self.Polyg[ip, :]))

            f.write('CELL_DATA %d\n' % num_dataset)
            f.write('SCALARS cell_scalars double\n')
            f.write('LOOKUP_TABLE default\n')
            for i in range(self.n_poly):
                f.write('%.15e\n' % strain[i])

        return

    def create_3Dmodel(self, theta_M: list | np.ndarray,
                       save_zip: bool = False):

        niter = len(theta_M)

        print('Create VTK...')
        # Delete previous data
        if os.path.exists(self.dir_save):
            shutil.rmtree(self.dir_save)
        os.makedirs(self.dir_save)
        # Generate vtk file
        from .vecmat.op_vecmat import calc_facet_area
        for ii in tqdm(range(niter)):
            vert_xyz = self.find_geo_miura_unit(theta_M[ii])
            # Dummy array for panel color
            farea = calc_facet_area(vert_xyz=vert_xyz, Polyg=self.Polyg)
            self.write_vtk(ii, vert_xyz, farea)

        if save_zip:
            zp = zipfile.ZipFile('%s.zip' % self.dir_save, 'w')
            dfile = glob.glob('%s/*.vtk' % self.dir_save)
            dfile = np.sort(dfile)
            for i in range(len(dfile)):
                zp.write(filename=dfile[i], arcname=None, compress_type=None, compresslevel=9)
            zp.close()

        return


if __name__ == "__main__":
    exit()
