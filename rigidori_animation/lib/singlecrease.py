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


class SingleCreaseAnalysis:
    def __init__(self, dir_save, fname_vtk):
        print('{0:*^49}'.format(''))
        print('* {0:^45} *'.format('Single Crease Static Analysis Tool'))
        print('{0:*^49}'.format(''))

        # VTK export
        self.dir_save = dir_save
        self.fname_vtk = fname_vtk

        pass

    # ******************************************************************
    #  Methods specific to single crease fold example
    # ******************************************************************

    def geo_init_simplefold(self, La: float, theta_M0: float, fig_out: bool = True):
        # **********************************************
        # Simple crease fold with two regular triangles
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
        self.n_node = 4
        self.n_edge = 5
        self.n_poly = 2

        # Define node
        vert_xyz = self.find_geo_simplefold(theta_M0, La)

        # Define creases
        EdgeConct = np.zeros((self.n_edge, 2), dtype=int)
        EdgeConct[0, :] = [0, 2]
        EdgeConct[1, :] = [0, 3]
        EdgeConct[2, :] = [1, 2]
        EdgeConct[3, :] = [1, 3]
        EdgeConct[4, :] = [2, 3]

        # Define polygon
        Polyg = np.zeros((self.n_poly, 3), dtype=int)
        Polyg[0, :] = [0, 2, 3]
        Polyg[1, :] = [1, 2, 3]

        self.vert_xyz0 = vert_xyz
        self.EdgeConct = EdgeConct
        self.Polyg = Polyg

        print('{0:*^49}'.format(''))
        print(' Initialize Simple Crease Fold Geometry')
        print('  {0:<24s} : {1:.6f}'.format('Side length', self.La))
        print('  {0:<24s} : {1:d} {2:d} {3:d}'.format('Node, crease, polygon', self.n_node, self.n_edge, self.n_poly))
        print('  {0:<24s} : {1:.6f} (deg)'.format('Initial fold angle', np.degrees(self.theta_M0)))
        print('{0:*^49}'.format(''))

        if fig_out:
            from .plot_geometry import plot_projection
            plot_projection(vert_xyz=vert_xyz,
                            Polyg=self.Polyg,
                            EdgeConct=self.EdgeConct,
                            figname='(SingleCrease-init)')

        return

    def find_geo_simplefold(self, theta_M: float, La: float):
        '''
        find_geo_simplefold
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
            f.write('POINTS %d double\n' % num_points)
            # Write x,y,z coodinates of all nodes
            for iv in range(num_points):
                f.write("%.15e %.15e %.15e\n" % (vert_xyz[iv, 0], vert_xyz[iv, 1], vert_xyz[iv, 2]))

            #   POLYGONS
            num_dataset = self.n_poly
            num_datanum = 4 * num_dataset
            f.write('POLYGONS %d %d\n' % (num_dataset, num_datanum))
            for ip in range(self.n_poly):
                f.write('3 %d %d %d\n' % (self.Polyg[ip, 0], self.Polyg[ip, 1], self.Polyg[ip, 2]))

            f.write('CELL_DATA %d\n' % num_dataset)
            f.write('SCALARS cell_scalars double\n')
            f.write('LOOKUP_TABLE default\n')
            for i in range(self.n_poly):
                f.write('%.15e\n' % strain[i])

        return

    def create_3Dmodel_simplefold(self, theta_M: npt.ArrayLike,
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
            vert_xyz = self.find_geo_simplefold(theta_M[ii], self.La)
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
