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
                       La: float, Lb: float, alpha: float,
                       theta_M0: float,
                       fig_out: bool = True):
        # **********************************************
        #        Miura-ori unit node definition
        # **********************************************
        #                    3
        #                  / | \
        #                 /  |  \
        #                /   |   \
        #               /    |    \
        #              /     |     \
        #             /      |      \
        #            /       |       \
        #           1        |        0
        #            \       |       /
        #             \      |      /
        #              \     |     /
        #               \    |    /
        #                \   |   /
        #                 \  |  /
        #                  \ | /
        #                    2
        # **********************************************
        self.La = La
        self.theta_M0 = theta_M0
        self.n_node = (2 * n_unitx + 3) * n_unity
        self.n_edge = 2 * n_unitx * (n_unity + 1) + 2 * n_unity * (n_unitx)
        self.n_poly = 4 * n_unitx * n_unity

        # Define nodal coordinates
        vert_xyz = self.find_geo_miura(theta_M0, La)

        # Define edge connection
        EdgeConct = np.zeros((self.n_edge, 2), dtype=int)

        # Define polygon
        Polyg = np.zeros((self.n_poly, 3), dtype=int)

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

    def find_geo_miura(self, theta_M: float, La: float):
        '''
        find_geo_miura
            Find the xyz coordinates of all nodes at specified theta_M.

        Args:
            theta_M (float): folding angle
            La (float): side length of regular triangle

        Returns:
            vert_xyz
        '''

        vert_xyz = np.zeros((self.n_node, 3))
        vert_xyz[0, 0] = 0.5 * np.sqrt(3) * La * np.cos(0.5 * (np.pi - theta_M))
        vert_xyz[0, 2] = 0.5 * np.sqrt(3) * La * np.sin(0.5 * (np.pi - theta_M))
        vert_xyz[1, 0] = -0.5 * np.sqrt(3) * La * np.cos(0.5 * (np.pi - theta_M))
        vert_xyz[1, 2] = 0.5 * np.sqrt(3) * La * np.sin(0.5 * (np.pi - theta_M))
        vert_xyz[2, 1] = -0.5 * La
        vert_xyz[3, 1] = 0.5 * La

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
            num_datanum = 4 * num_dataset
            f.write('POLYGONS %d %d\n' % (num_dataset, num_datanum))
            for ip in range(self.n_poly):
                f.write('3 %d %d %d\n' % (self.Polyg[ip, 0], self.Polyg[ip, 1], self.Polyg[ip, 2]))

            f.write('CELL_DATA %d\n' % num_dataset)
            f.write('SCALARS cell_scalars float\n')
            f.write('LOOKUP_TABLE default\n')
            for i in range(self.n_poly):
                f.write('%e\n' % strain[i])

        return

    def create_3Dmodel_simplefold(self, theta_M: list | np.ndarray):

        niter = len(theta_M)
        # VTK export
        self.dir_save = './vtk_singlecrease'
        self.fname_vtk = 'singlecrease'

        print('Create VTK...')
        # Delete previous data
        if os.path.exists(self.dir_save):
            shutil.rmtree(self.dir_save)
        os.makedirs(self.dir_save)
        # Generate vtk file
        for ii in tqdm(range(niter)):
            vert_xyz = self.find_geo_simplefold(theta_M[ii], self.La)
            # Dummy array for panel color
            strain = (theta_M[ii] - theta_M[0]) / (theta_M[-1] - theta_M[0]) * np.ones(self.n_poly)
            self.write_vtk(ii, vert_xyz, strain)

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
