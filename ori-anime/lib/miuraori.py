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

        self.La = La
        self.Lb = Lb
        self.Ld = Ld
        self.alpha = alpha
        self.theta_M0 = theta_M0
        n_node = (n_unitx + 1) * (n_unity + 1)
        n_poly = n_unitx * n_unity
        n_elem = n_unitx * (n_unity + 1) + n_unity * (n_unitx + 1)

        print('{0:*^49}'.format(''))
        print(' Initialize Miura-ori Fold Geometry')
        print('  {0:<24s} : {1:.6f}'.format('Side length (a)', self.La))
        print('  {0:<24s} : {1:.6f}'.format('Side length (b)', self.Lb))
        print('  {0:<24s} : {1:.6f}'.format('Apex angle (alpha)', np.degrees(self.alpha)))
        print('  {0:<24s} : {1:d} {2:d}'.format('Nx, Ny', n_unitx, n_unity))
        print('  {0:<24s} : {1:d} {2:d} {3:d}'.format('Node, crease, polygon', n_node, n_elem, n_poly))
        print('  {0:<24s} : {1:.6f} (deg)'.format('Initial fold angle', np.degrees(self.theta_M0)))
        print('{0:*^49}'.format(''))

        # Define polygon
        Polyg = np.zeros((n_poly, 4), dtype=int)
        for iy in range(n_unity):
            for ix in range(n_unitx):
                Polyg[ix + iy * n_unitx, :] = [iy * (n_unitx + 1) + ix, iy * (n_unitx + 1) + ix + 1,
                                               (iy + 1) * (n_unitx + 1) + ix + 1, (iy + 1) * (n_unitx + 1) + ix]
        # Define edge connection
        Elem = np.zeros((n_elem, 2), dtype=int)
        ii = 0
        for iy in range(n_unity + 1):
            for ix in range(n_unitx):
                Elem[ii, :] = [iy * (n_unitx + 1) + ix, iy * (n_unitx + 1) + ix + 1]
                ii = ii + 1
        for ix in range(n_unitx + 1):
            for iy in range(n_unity):
                Elem[ii, :] = [(n_unitx + 1) * iy + ix, (n_unitx + 1) * (iy + 1) + ix]
                ii = ii + 1

        # Define nodal coordinates
        vert_xyz = self.find_geo_miura(theta_M0,
                                       n_node, n_unitx, n_unity)

        if fig_out:
            from .plot_geometry import plot_projection
            plot_projection(vert_xyz=vert_xyz,
                            Polyg=Polyg,
                            EdgeConct=Elem,
                            figname='(Miura-init)')

        self.vert_xyz0 = vert_xyz
        self.Elem = Elem
        self.Polyg = Polyg
        self.n_unitx = n_unitx
        self.n_unity = n_unity
        self.n_node = n_node
        self.n_elem = n_elem
        self.n_poly = n_poly

        return

    def geo_vec(self, theta_M):
        vec_base = {}

        ct = np.cos(theta_M) * np.tan(self.alpha)
        ss = np.sin(theta_M) * np.sin(self.alpha)
        vec_base['e0'] = np.array([0.,
                                   self.Lb * np.sqrt(1 - ss**2),
                                   0.])
        vec_base['e1'] = (self.La / np.sqrt(1 + ct**2)) * np.array([0.,
                                                                    1.,
                                                                    ct])
        vec_base['e1p'] = (self.La / np.sqrt(1 + ct**2)) * np.array([0.,
                                                                     1.,
                                                                     -ct])
        vec_base['e2'] = self.Lb * np.array([ss,
                                             -np.sqrt(1. - ss**2),
                                             0.])
        vec_base['e2p'] = self.Lb * np.array([ss,
                                             np.sqrt(1. - ss**2),
                                             0.])

        return vec_base

    def find_geo_miura(self, theta_M: float,
                       n_node: int,
                       n_unitx: int = 1, n_unity: int = 1):
        ex = np.array([1, 0, 0])
        ey = np.array([0, 1, 0])
        ez = np.array([0, 0, 1])
        vert_xyz = np.zeros((n_node, 3))
        vec_base = self.geo_vec(theta_M)
        e0 = vec_base['e0']
        e1 = vec_base['e1']
        e2 = vec_base['e2']
        e1p = vec_base['e1p']
        e2p = vec_base['e2p']
        ii = 0

        vert_xyz[ii] = [0, 0, 0]
        ii = ii + 1
        for ix in range(1, n_unitx + 1):
            if np.mod(ix, 2) == 1:
                vert_xyz[ii] = vert_xyz[ii - 1] + self.Lb * e2p
                ii = ii + 1
            elif np.mod(ix, 2) == 0:
                vert_xyz[ii] = vert_xyz[ii - 1] + self.Lb * e2
                ii = ii + 1
        for iy in range(1, n_unity + 1):
            if np.mod(iy, 2) == 1:
                # print(ii, ii + (n_unitx + 1))
                vert_xyz[ii:ii + (n_unitx + 1)] = vert_xyz[ii - (n_unitx + 1):ii] + self.La * e1
                ii = ii + (n_unitx + 1)
            elif np.mod(iy, 2) == 0:
                vert_xyz[ii:ii + (n_unitx + 1)] = vert_xyz[ii - (n_unitx + 1):ii] + self.La * e1p
                ii = ii + (n_unitx + 1)

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
            # vert_xyz = self.find_geo_miura_unit(theta_M[ii])
            vert_xyz = self.find_geo_miura(theta_M[ii],
                                           self.n_node,
                                           self.n_unitx, self.n_unity)
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
