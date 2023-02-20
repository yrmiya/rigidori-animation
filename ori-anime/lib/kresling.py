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


class KreslingAnalysis:
    def __init__(self, dir_save, fname_vtk) -> None:
        self.dir_save = dir_save
        self.fname_vtk = fname_vtk
        pass

    def calc_creaselength(self, h0, radius, angle0, n_sec):
        a0 = np.sqrt(h0**2 + 4. * radius**2 * np.sin(0.5 * angle0 - 0.5 * np.pi / n_sec)**2)
        b0 = np.sqrt(h0**2 + 4. * radius**2 * np.sin(0.5 * angle0 + 0.5 * np.pi / n_sec)**2)

        return a0, b0

    def geo_init_kres(self, n_unit, n_side, Lc, h0, ph0, fig_out=False):
        # *********************************************************************
        # *********************************************************************
        #  (N-2)-[N_i-1]-(N-1)
        #    |           / |
        #    |         /   |
        #    |       /     |
        #    |   [N_i-2]   |
        #    |     /       |
        #    |   /         |
        #    | /           |
        #  (N-4)-[N_i-3]-(N-3)
        #    |           / |
        #    :      :      :
        #    | /           |
        #  ( n )---[2]---(n+1)---[0]---(n+2)---...---(2n-1)---[-]---( n )
        #    |           / |
        #    |         /   |
        #    |       /     |
        #  [N_i]   [1]     |
        #    |     /       |
        #    |   /         |
        #    | /           |
        #  ( 0 )---[0]---( 1 )---[0]---( 2 )---...---( n-1)---[-]---( 0 )
        #
        # *********************************************************************
        # *********************************************************************

        n_node = (n_unit + 1) * 2
        n_poly = (n_unit) * 2
        th_c = 2.0 * np.pi / n_side
        Lr = 0.5 * Lc / np.sin(0.5 * th_c)

        n_node = n_side * (n_unit + 1)  # (main)
        n_poly = n_side * (n_unit) + (n_unit + 1)  # (main) + (caps)
        n_poly_main = n_side * (n_unit)   # (main) + (caps)
        n_poly_cap = n_unit + 1  # (main) + (caps)

        # ********************
        #      Vertices
        # ********************
        vert_xyz = np.zeros((n_node, 3))

        for i in range(n_unit + 1):
            print(n_side * i, n_side * (i + 1))
            vert_xyz[n_side * i:n_side * (i + 1), 2] = h0 * i
            for j in range(n_side):
                vert_xyz[n_side * i + (j - 1), 0] = Lr * np.cos(-0.5 * (np.pi + th_c) + j * th_c + i * (ph0 - 0.5 * th_c))
                vert_xyz[n_side * i + (j - 1), 1] = Lr * np.sin(-0.5 * (np.pi + th_c) + j * th_c + i * (ph0 - 0.5 * th_c))

        # ********************
        #      Polygons
        # ********************
        # Constituent vertices
        Polyg = np.zeros((n_poly), dtype=object)
        for i in range(n_unit):
            for j in range(n_side):
                n0 = j + n_side * i
                if j == 0:
                    Polyg[2 * n_side * i + 2 * j] = [n0, n0 + n_side, n0 + 2 * n_side - 1]
                    Polyg[2 * n_side * i + 2 * j + 1] = [n0, n0 + 1, n0 + n_side]
                elif j == n_side - 1:
                    Polyg[2 * n_side * i + 2 * j] = [n0, n0 + n_side, n0 + n_side - 1]
                    Polyg[2 * n_side * i + 2 * j + 1] = [n0, n0 - n_side + 1, n0 + n_side]
                else:
                    Polyg[2 * n_side * i + 2 * j] = [n0, n0 + n_side, n0 + n_side - 1]
                    Polyg[2 * n_side * i + 2 * j + 1] = [n0, n0 + 1, n0 + n_side]
        for i in range(n_unit + 1):
            Polyg[-(i + 1)] = np.arange(n_side * i, n_side * (i + 1), 1)

        # from .plot_geometry import plot_projection
        # plot_projection(vert_xyz=vert_xyz,
        #                 Polyg=Polyg,
        #                 EdgeConct=EdgeConct,
        #                 figname='(Miura-init)')

        return vert_xyz, Polyg

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
            f.write('Kresling\n')
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

if __name__ == "__main__":
    h0 = 30e-3
    radius = 36e-3
    angle0 = np.radians(70)
    n_sec = 6
    koa = KreslingAnalysis()
    a0, b0 = koa.calc_creaselength(h0, radius, angle0, n_sec)
    print(a0, b0)
    exit()
