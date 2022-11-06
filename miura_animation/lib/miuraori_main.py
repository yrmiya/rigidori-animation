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
plt.style.use('./common/custom.mplstyle')


class MiuraOriAnalysis:
    def __init__(self,):
        print('{0:*^49}'.format(''))
        print('* {0:^45} *'.format('Miura-ori Static Analysis Tool'))
        print('{0:*^49}'.format(''))

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
            self.plot_projection(vert_xyz, figname='(init)')

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

    def calc_facet_area(self, vert_xyz: npt.ArrayLike):
        facet_area = np.zeros(self.n_poly)
        for ip in range(self.n_poly):
            vec1 = vert_xyz[self.Polyg[ip, 0], :] - vert_xyz[self.Polyg[ip, 1], :]
            vec2 = vert_xyz[self.Polyg[ip, 2], :] - vert_xyz[self.Polyg[ip, 1], :]
            facet_area[ip] = np.linalg.norm(np.cross(vec1, vec2))

        return facet_area

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
            f.write('Single crease fold\n')
            f.write('ASCII\n')
            num_points = np.size(vert_xyz, axis=0)
            f.write('DATASET POLYDATA\n')
            # Nodes
            # Specify number of nodes
            f.write('POINTS %d float\n' % num_points)
            # Write x,y,z coodinates of all nodes
            for iv in range(num_points):
                f.write("%f %f %f\n" % (vert_xyz[iv, 0], vert_xyz[iv, 1], vert_xyz[iv, 2]))

            #   POLYGONS
            num_dataset = self.n_poly
            num_datanum = 5 * num_dataset
            f.write('POLYGONS %d %d\n' % (num_dataset, num_datanum))
            for ip in range(self.n_poly):
                f.write('4 %d %d %d %d\n' % tuple(self.Polyg[ip, :]))

            f.write('CELL_DATA %d\n' % num_dataset)
            f.write('SCALARS cell_scalars float\n')
            f.write('LOOKUP_TABLE default\n')
            for i in range(self.n_poly):
                f.write('%e\n' % strain[i])

        return

    def create_3Dmodel(self, theta_M: list | np.ndarray,
                       save_zip: bool = False):

        niter = len(theta_M)
        # VTK export
        self.dir_save = './vtk_miura'
        self.fname_vtk = 'miura'

        print('Create VTK...')
        # Delete previous data
        if os.path.exists(self.dir_save):
            shutil.rmtree(self.dir_save)
        os.makedirs(self.dir_save)
        # Generate vtk file
        for ii in tqdm(range(niter)):
            vert_xyz = self.find_geo_miura_unit(theta_M[ii])
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

    def plot_projection(self, vert_xyz: list | np.ndarray, figname: str = ''):
        '''
        plot_projection
            Plot xy-plane, xz-plane, and yz-plane projection of the origami

        Args:
            vert_xyz (list | np.ndarray): xyz coordinates
            figname (str, optional): Name of the figure. Defaults to ''.
        '''

        fig, ax = plt.subplots(1, 3, num='xy projection %s' % (figname), figsize=(18, 6))
        # xy plane
        # for i in range(self.n_node):
        for i in range(self.n_node):
            ax[0].scatter(vert_xyz[i, 0], vert_xyz[i, 1], s=80, zorder=2,)
        for i in range(self.n_edge):
            ax[0].plot(vert_xyz[self.EdgeConct[i, :], 0],
                       vert_xyz[self.EdgeConct[i, :], 1],
                       color='#444444',
                       linewidth=2,
                       zorder=1,
                       )
        for i in range(self.n_poly):
            ax[0].fill(vert_xyz[self.Polyg[i, :], 0],
                       vert_xyz[self.Polyg[i, :], 1],
                       #    color='#444444',
                       alpha=0.2,
                       linewidth=2,
                       zorder=0,
                       )
        ax[0].set_xlabel('$x$')
        ax[0].set_ylabel('$y$')
        ax[0].set_aspect('equal', 'box')
        # xz plane
        for i in range(self.n_node):
            ax[1].scatter(vert_xyz[i, 0], vert_xyz[i, 2], s=80, zorder=2,)
        for i in range(self.n_edge):
            ax[1].plot(vert_xyz[self.EdgeConct[i, :], 0],
                       vert_xyz[self.EdgeConct[i, :], 2],
                       color='#444444',
                       linewidth=2,
                       zorder=1,
                       )
        for i in range(self.n_poly):
            ax[1].fill(vert_xyz[self.Polyg[i, :], 0],
                       vert_xyz[self.Polyg[i, :], 2],
                       #    color='#444444',
                       alpha=0.2,
                       linewidth=2,
                       zorder=0,
                       )
        ax[1].set_xlabel('$x$')
        ax[1].set_ylabel('$z$')
        ax[1].set_aspect('equal', 'box')
        # yz plane
        for i in range(self.n_node):
            ax[2].scatter(vert_xyz[i, 1], vert_xyz[i, 2], s=80, zorder=2,)
        for i in range(self.n_edge):
            ax[2].plot(vert_xyz[self.EdgeConct[i, :], 1],
                       vert_xyz[self.EdgeConct[i, :], 2],
                       color='#444444',
                       linewidth=2,
                       zorder=1,
                       )
        for i in range(self.n_poly):
            ax[2].fill(vert_xyz[self.Polyg[i, :], 1],
                       vert_xyz[self.Polyg[i, :], 2],
                       #    color='#444444',
                       alpha=0.2,
                       linewidth=2,
                       zorder=0,
                       )
        ax[2].set_xlabel('$y$')
        ax[2].set_ylabel('$z$')
        ax[2].set_aspect('equal', 'box')
        fig.tight_layout()

        return


if __name__ == "__main__":
    exit()
