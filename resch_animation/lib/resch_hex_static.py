# -*- coding: utf-8 -*-
"""
Written by: Yasuhiro Miyazawa
"""

# Import packages
import numpy as np
import numpy.typing as npt
from scipy import optimize
from scipy import interpolate
import matplotlib
matplotlib.use('QtAgg')  # Backend
import matplotlib.pyplot as plt

import sys
import os
import shutil
import zipfile

from tqdm import tqdm

# Figure parameters
plt.style.use('./common/custom.mplstyle')


class ReschOrigamiAnalysis:
    def __init__(self,):
        print('{0:*^49}'.format(''))
        print('* {0:^45} *'.format('Resch Origami Static Analysis Tool'))
        print('{0:*^49}'.format(''))

        # axis 60 rotation about z on hexagon
        theta1 = np.pi / 3
        T_4 = np.array([[np.cos(theta1), -np.sin(theta1), 0], [np.sin(theta1), np.cos(theta1), 0], [0, 0, 1]])
        # axis 120 rotation about z on hexagon
        theta2 = 2 * np.pi / 3
        T_5 = np.array([[np.cos(theta2), np.sin(theta2), 0], [-np.sin(theta2), np.cos(theta2), 0], [0, 0, 1]])
        T_6 = T_4.T

        self.T_4 = T_4
        self.T_5 = T_5
        self.T_6 = T_6

        pass

    def geo_init(self, La: float, n_orbit: int, theta_M0: float):
        n_hex = 1 + 3 * (n_orbit + n_orbit**2)
        n_node_hex = int(6 * n_hex)
        n_node_hex_origin = int(n_hex)
        n_node_tri = int(6 * (n_orbit + n_orbit**2))
        n_node = n_node_hex + n_node_tri + n_node_hex_origin
        n_edge_hex = int(6 * n_hex)
        n_edge_tri = int(7 * n_node_tri)
        n_edge = n_edge_hex + n_edge_tri
        n_poly_hex = n_hex
        n_poly_tri = int(6 * (n_node_tri - n_orbit * 3))
        n_poly = n_poly_hex + n_poly_tri

        # Determine crease length
        Lb = La / np.sqrt(3.)
        Lc = 2. * La / np.sqrt(3.)

        self.La = La
        self.Lb = Lb
        self.Lc = Lc
        self.n_orbit = n_orbit
        self.theta_M0 = theta_M0
        #
        self.n_node_hex = n_node_hex
        self.n_node_tri = n_node_tri
        self.n_node = n_node
        #
        self.n_edge_hex = n_edge_hex
        self.n_edge_tri = n_edge_tri
        self.n_edge = n_edge
        #
        self.n_poly_hex = n_poly_hex
        self.n_poly_tri = n_poly_tri
        self.n_poly = n_poly

        # Import fold angles
        data_th = np.genfromtxt('./lookuptable/theta_MKJS.csv', delimiter=',')
        theta_M = data_th[:, 0]
        theta_K = data_th[:, 1]
        theta_J = data_th[:, 2]
        theta_S = data_th[:, 3]
        f_MK = interpolate.interp1d(theta_M, theta_K, kind='cubic')
        f_MJ = interpolate.interp1d(theta_M, theta_J, kind='cubic')
        f_MS = interpolate.interp1d(theta_M, theta_S, kind='cubic')

        theta_K0 = f_MK(theta_M0)
        theta_J0 = f_MJ(theta_M0)
        theta_S0 = f_MS(theta_M0)

        if n_orbit == 1:
            vert_xyz, EdgeConct, Polyg_tri, Polyg_hex,\
                PolygAdj, Polyg2Edge = self.geo_init_1orbit(La, theta_M0)

        print('{0:*^49}'.format(''))
        print(' Initialize Resch tessellation')
        print('  {0:<32s} : {1:<12d}'.format('Number of orbit', n_orbit))
        print('  {0:<32s} : {1:<.12f}'.format('Hexagon side length', La))
        print('  {0:<32s} : {1:<.12f}'.format('Minor crease length', Lb))
        print('  {0:<32s} : {1:<.12f}'.format('Major crease length', Lc))
        print('  {0:<32s}'.format('Initial fold angles th_M0'))
        print('    {0:<30s} : {1:<.8f} (deg)'.format('th_M0', np.degrees(theta_M0)))
        print('    {0:<30s} : {1:<.8f} (deg)'.format('th_K0', np.degrees(theta_K0)))
        print('    {0:<30s} : {1:<.8f} (deg)'.format('th_J0', np.degrees(theta_J0)))
        print('    {0:<30s} : {1:<.8f} (deg)'.format('th_S0', np.degrees(theta_S0)))
        print('  {0:<32s}'.format('Node, crease, polygon'))
        print('    {0:<30s} : {1:d}, {2:d}, {3:d}'.format('Total', n_node, n_edge, n_poly))
        print('      {0:<28s} : {1:d}, {2:d}, {3:d}'.format('Hexagon', n_node_hex, n_edge_hex, n_poly_hex))
        print('      {0:<28s} : {1:d}, {2:d}, {3:d}'.format('Triangle', n_node_tri, n_edge_tri, n_poly_tri))
        print('{0:*^49}'.format(''))

        return

    def geo_init_1orbit(self, La, theta_M0):

        self.La = La
        # self.n_node = 54

        # n_tri = 48  # 6*6+2*6=48
        # n_hex = 7
        n_tri = self.n_poly_tri
        n_hex = self.n_poly_hex

        EdgeConct = np.zeros((self.n_edge, 2), dtype=int)
        PolygAdj = np.zeros((self.n_edge, 2), dtype=int)
        Polyg2Edge = np.zeros((self.n_poly, 2), dtype=int)

        # Define triangles
        Polyg_tri = np.zeros((n_tri, 3), dtype=int)

        Polyg_tri[0, :] = [0, 1, 43]
        Polyg_tri[1, :] = [1, 2, 44]
        Polyg_tri[2, :] = [2, 3, 45]
        Polyg_tri[3, :] = [3, 4, 46]
        Polyg_tri[4, :] = [4, 5, 47]
        Polyg_tri[5, :] = [5, 0, 42]
        Polyg_tri[6, :] = [0, 9, 42]
        Polyg_tri[7, :] = [1, 16, 43]
        Polyg_tri[8, :] = [2, 23, 44]
        Polyg_tri[9, :] = [3, 24, 45]
        Polyg_tri[10, :] = [4, 31, 46]
        Polyg_tri[11, :] = [5, 38, 47]
        Polyg_tri[12, :] = [0, 9, 43]
        Polyg_tri[13, :] = [1, 16, 44]
        Polyg_tri[14, :] = [2, 23, 45]
        Polyg_tri[15, :] = [3, 24, 46]
        Polyg_tri[16, :] = [4, 31, 47]
        Polyg_tri[17, :] = [5, 38, 42]
        Polyg_tri[18, :] = [8, 9, 43]
        Polyg_tri[19, :] = [15, 16, 44]
        Polyg_tri[20, :] = [22, 23, 45]
        Polyg_tri[21, :] = [24, 29, 46]
        Polyg_tri[22, :] = [30, 31, 47]
        Polyg_tri[23, :] = [37, 38, 42]
        Polyg_tri[24, :] = [8, 17, 43]
        Polyg_tri[25, :] = [16, 17, 43]
        Polyg_tri[26, :] = [15, 18, 44]
        Polyg_tri[27, :] = [18, 23, 44]
        Polyg_tri[28, :] = [22, 25, 45]
        Polyg_tri[29, :] = [24, 25, 45]
        Polyg_tri[30, :] = [29, 32, 46]
        Polyg_tri[31, :] = [31, 32, 46]
        Polyg_tri[32, :] = [30, 39, 47]
        Polyg_tri[33, :] = [38, 39, 47]
        Polyg_tri[34, :] = [10, 37, 42]
        Polyg_tri[35, :] = [9, 10, 42]
        Polyg_tri[36, :] = [8, 17, 48]
        Polyg_tri[37, :] = [12, 17, 48]
        Polyg_tri[38, :] = [15, 18, 49]
        Polyg_tri[39, :] = [18, 19, 49]
        Polyg_tri[40, :] = [25, 26, 50]
        Polyg_tri[41, :] = [22, 25, 50]
        Polyg_tri[42, :] = [32, 33, 51]
        Polyg_tri[43, :] = [29, 32, 51]
        Polyg_tri[44, :] = [30, 39, 52]
        Polyg_tri[45, :] = [39, 40, 52]
        Polyg_tri[46, :] = [10, 11, 53]
        Polyg_tri[47, :] = [10, 37, 53]

        # Define hexagons
        Polyg_hex = np.zeros((n_hex, 6), dtype=int)
        for i in range(self.n_poly_hex):
            Polyg_hex[i, :] = np.arange(6 * i, 6 * (i + 1), dtype=int)

        self.Polyg_tri = Polyg_tri
        self.Polyg_hex = Polyg_hex
        self.n_poly_tri = n_tri
        self.n_poly_hex = n_hex

        theta_M = np.array([theta_M0, ])
        vert_xyz = self.solve_geo_1orbit(theta_M=theta_M, vtk_out=False, fig_out=False)

        self.Polyg = np.zeros((self.n_poly, 6), dtype=int)
        self.Polyg[:self.n_poly_hex, :] = Polyg_hex
        # self.Polyg[self.n_poly_hex:, :3] = Polyg_tri

        return vert_xyz, EdgeConct, Polyg_tri, Polyg_hex, PolygAdj, Polyg2Edge

    def calc_ukvector(self, theta_K: float, theta_M: float, La: float):
        u1_b = La * np.array([0.5 * np.sin(theta_M), 0.5 / np.sqrt(3.), 0.5 * np.cos(theta_M)])
        u2_b = La * np.array([-0.5 * np.sin(theta_M), 0.5 * np.sqrt(3.), -0.5 * np.cos(theta_M)])

        n1 = u2_b / La

        mu1 = 2. * theta_K
        k1_b = np.dot(u1_b, n1) * n1 * (1. - np.cos(mu1)) + np.cos(mu1) * u1_b - np.sin(mu1) * np.cross(n1, u1_b)

        E = u1_b + u2_b - k1_b
        G = u1_b
        EG = G - E
        # n2 = EG / (2. * a / np.sqrt(3.))
        n2 = 0.5 * EG * np.sqrt(3) / La
        mu2 = 2 * theta_M
        k4_b = np.dot(k1_b, n2) * n2 * (1 - np.cos(mu2)) + np.cos(mu2) * k1_b - np.sin(mu2) * np.cross(n2, k1_b)
        k3_b = -(u2_b - k1_b + k4_b)
        u3_b = La * np.array([0.5 * np.sin(theta_M), 0.5 * np.sqrt(3.), -0.5 * np.cos(theta_M)])

        psi = np.pi - np.arccos(np.sin(theta_M)) - np.arccos(1. / (6. * np.sqrt(1. / 9. - 1. / 12. * np.sin(theta_M)**2))) - \
            np.emath.arccos((2. / 3. - np.sin(theta_M)**2) / (2. * np.sqrt(1. / 9. - 1. / 12. * np.sin(theta_M)**2))).real

        T_1 = np.array([[1., 0., 0.], [0., np.cos(psi), -np.sin(psi)], [0., np.sin(psi), np.cos(psi)]])
        t_1 = T_1.T
        p1 = 2. * np.pi / 3.
        T_3 = np.array([[np.cos(p1), np.sin(p1), 0.], [-np.sin(p1), np.cos(p1), 0.], [0., 0., 1.]])
        u4_b = np.dot(np.dot(np.dot(t_1, T_3), T_1), u3_b)

        return psi, u1_b, u2_b, u3_b, u4_b, k1_b, k3_b, k4_b

    def funcMK(self, theta_K: float, theta_M: float, La: float):

        _, _, _, _, u4_b, _, k3_b, _ = self.calc_ukvector(theta_K, theta_M, La)

        f = np.dot(u4_b, k3_b) - La**2 * np.cos(np.pi / 3.)

        return f

    def calc_theta_JS(self, theta_M: float, theta_K: float, La: float):

        psi, u1_b, _, _, u4_b, _, k3_b, _ = self.calc_ukvector(theta_K, theta_M, La)

        v1_b = k3_b - 0.5 * u4_b
        theta_J = np.emath.arccos(2. * np.dot(u1_b, v1_b) / La**2).real
        theta_S = np.arcsin(np.cos(psi))

        return theta_J, theta_S

    def find_theta_MKJS(self, theta_M: list | np.ndarray,
                        La: float = 1.0,
                        x0: float = 0.0,
                        method: str = 'hybr',
                        fig_out: bool = True,
                        save_lookuptab: bool = False
                        ):
        '''
        find_theta_MKJS
        Determines theta_K and theta_J as a function of theta_M

        Args:
            theta_M (list | np.ndarray): Array of theta_M
            a (float, optional): Side length. Defaults to 1.0.
            x0 (float, optional): Initial guess for the first step. Defaults to 0.0.
            method (str, optional): Method for scipy.optimize.root. Defaults to 'hybr'.
            fig_out (bool, optional): If True, displays figures. Defaults to True.
            save_lookuptab (bool, optional): If True, saves lookup table. Defaults to False.

        Returns:
            (list | np.ndarray): List of theta_K and theta_J
        '''

        niter = len(theta_M)

        # Allocate arrays to store solutions
        theta_J = np.zeros_like(theta_M)
        theta_K = np.zeros_like(theta_M)
        theta_S = np.zeros_like(theta_M)

        for i in tqdm(range(niter)):
            # Find root of funcMK at each theta_M
            output = optimize.root(self.funcMK, x0, args=(theta_M[i], La),
                                   method=method)
            # Save solution
            theta_K[i] = output.x
            # Update initial guess for next step
            x0 = theta_K[i]

        # Substitute first and last entry with correct mathematical values
        theta_K[0] = 0.0
        theta_K[-1] = 0.5 * np.pi

        # Calculate theta_J and theta_S
        for i in range(niter):
            theta_J[i], theta_S[i] = self.calc_theta_JS(theta_M=theta_M[i], theta_K=theta_K[i], La=La)

        # Store fold angles
        theta_MKJS = np.zeros((4, len(theta_M)))
        theta_MKJS[0, :] = theta_M
        theta_MKJS[1, :] = theta_K
        theta_MKJS[2, :] = theta_J
        theta_MKJS[3, :] = theta_S

        # Save lookup table
        if save_lookuptab:
            print(' Save lookup table')
            fpath = './lookuptable/theta_MKJS.csv'
            file_exists = os.path.exists(fpath)
            if file_exists:
                with open(fpath, 'r') as fp:
                    for count, line in enumerate(fp):
                        pass
                count = count + 1
                file_stat = os.stat(fpath)
                if file_stat.st_size >= 1e3 and file_stat.st_size < 1e6:
                    fsize = file_stat.st_size * 1e-3
                    funit = 'KB'
                elif file_stat.st_size >= 1e6:
                    fsize = file_stat.st_size * 1e-6
                    funit = 'MB'
                else:
                    fsize = file_stat.st_size
                    funit = 'B'
                print('  ! File already exists ({0:d} lines; {1:.2f} {2:2s}) !'.format(count, fsize, funit))
                value = input('  Overwrite ? [Press Y to overwrite, press elsewhere otherwise]\n')
                if value == 'Y' or value == 'y':
                    print('   Saving lookup table...')
                    data = np.stack((theta_M, theta_K, theta_J, theta_S), axis=-1)
                    np.savetxt(fpath, data, delimiter=',')
                    print('    ...Saved')
                else:
                    print('   File not saved.')
                    pass
            else:
                print('   Saving lookup table...')
                data = np.stack((theta_M, theta_K, theta_J, theta_S), axis=-1)
                np.savetxt(fpath, data, delimiter=',')
                print('    ...Saved')

        # Generate figures
        if fig_out:
            plt.figure('theta MKJS')
            plt.plot(np.degrees(theta_M), np.degrees(theta_K), label='$\\theta_K$')
            plt.plot(np.degrees(theta_M), np.degrees(theta_J), label='$\\theta_J$')
            plt.plot(np.degrees(theta_M), np.degrees(theta_S), label='$\\theta_S$')
            plt.axhline(y=90.0, linewidth=1.0, color='#484848', linestyle='dashed')
            plt.legend(loc=4)
            plt.xlim(np.degrees(min(theta_M)), np.degrees(max(theta_M)))
            plt.ylim(0, 180)
            plt.xlabel('$\\theta_M$ (deg)')
            plt.ylabel('$\\theta_K$, $\\theta_J$ (deg)')

        return theta_MKJS

    def BaseHexLayer(self, theta_M, theta_J):
        A = np.zeros((8, 6, 3))
        B = np.zeros((8, 6, 3))
        C = np.zeros((8, 6, 3))
        D = np.zeros((8, 6, 3))
        E = np.zeros((8, 6, 3))

        A[0, 1, :] = np.array([0, self.La, 0])
        A[0, 0, :] = np.array([-0.5 * np.sqrt(3) * self.La, 0.5 * self.La, 0])
        A[0, 5, :] = np.array([-0.5 * np.sqrt(3) * self.La, -0.5 * self.La, 0])

        A1B1 = self.La / np.sqrt(3) * np.array([1 / 2 * np.cos(theta_J), -np.sqrt(3) / 2 * np.cos(theta_J), -np.sin(theta_J)])
        B[0, 1, :] = A[0, 1, :] + A1B1
        B1A1 = -A1B1
        A0B1 = B[0, 1] - A[0, 0]
        n3 = A0B1 / np.linalg.norm(A0B1)
        mu3 = 2 * theta_M
        B1C1 = np.dot(B1A1, n3) * n3 * (1 - np.cos(mu3)) + np.cos(mu3) * B1A1 - np.sin(mu3) * np.cross(n3, B1A1)
        C[0, 1] = B[0, 1] + B1C1
        C1A0 = A[0, 0] - C[0, 1]
        n4 = B1C1 / np.linalg.norm(B1C1)
        psi = np.pi - np.arccos(np.sin(theta_M)) - np.arccos(1. / (6 * (1 / 9 - 1 / 12 * np.sin(theta_M)**2)**0.5)) - \
            np.emath.arccos((2 / 3 - np.sin(theta_M)**2) / (2 * (1 / 9 - 1 / 12 * np.sin(theta_M)**2)**0.5)).real
        thetaS = np.arcsin(np.cos(psi))
        mu4 = 2 * thetaS
        C1E1 = np.dot(C1A0, n4) * n4 * (1 - np.cos(mu4)) + np.cos(mu4) * C1A0 - np.sin(mu4) * np.cross(n4, C1A0)
        E[0, 1] = C[0, 1] + C1E1
        E1B1 = B[0, 1] - E[0, 1]
        n5 = E1B1 / np.linalg.norm(E1B1)
        mu5 = mu3
        B1D1 = np.dot(B1C1, n5) * n5 * (1 - np.cos(mu5)) + np.cos(mu5) * B1C1 - np.sin(mu5) * np.cross(n5, B1C1)
        D[0, 1] = B[0, 1] + B1D1

        B[0, 0] = self.T_4 @ B[0, 1]
        C[0, 0] = self.T_4 @ C[0, 1]
        E[0, 0] = self.T_4 @ E[0, 1]
        D[0, 0] = self.T_4 @ D[0, 1]
        A[0, 2] = self.T_5 @ A[0, 0]
        B[0, 2] = self.T_5 @ B[0, 0]
        C[0, 2] = self.T_5 @ C[0, 0]
        E[0, 2] = self.T_5 @ E[0, 0]
        D[0, 2] = self.T_5 @ D[0, 0]
        A[0, 3] = self.T_5 @ A[0, 1]
        B[0, 3] = self.T_5 @ B[0, 1]
        C[0, 3] = self.T_5 @ C[0, 1]
        E[0, 3] = self.T_5 @ E[0, 1]
        D[0, 3] = self.T_5 @ D[0, 1]
        A[0, 4] = self.T_5 @ A[0, 2]
        B[0, 4] = self.T_5 @ B[0, 2]
        C[0, 4] = self.T_5 @ C[0, 2]
        E[0, 4] = self.T_5 @ E[0, 2]
        D[0, 4] = self.T_5 @ D[0, 2]
        A[0, 5] = self.T_5 @ A[0, 3]
        B[0, 5] = self.T_5 @ B[0, 3]
        C[0, 5] = self.T_5 @ C[0, 3]
        E[0, 5] = self.T_5 @ E[0, 3]
        D[0, 5] = self.T_5 @ D[0, 3]

        return A, B, C, D, E

    def M01_1(self, angles, A, C, D, E):
        r1 = np.array([[1, 0, 0], [0, np.cos(angles[0]), np.sin(angles[0])], [0, -np.sin(angles[0]), np.cos(angles[0])]])
        r2 = np.array([[np.cos(angles[1]), 0, -np.sin(angles[1])], [0, 1, 0], [np.sin(angles[1]), 0, np.cos(angles[1])]])
        r3 = r2 @ r1

        N1 = np.cross((A[0, 2] - A[0, 3]), (A[0, 4] - A[0, 3]))
        N2 = np.cross((E[0, 1] - C[0, 1]), (D[0, 0] - C[0, 1]))
        f = r3 @ N1
        f = f - N2

        return f[0:2]

    def M01_2(self, angles, Rangles, A, C, E):
        R1 = np.array([[1, 0, 0], [0, np.cos(Rangles[0]), np.sin(Rangles[0])], [0, -np.sin(Rangles[0]), np.cos(Rangles[0])]])
        R2 = np.array([[np.cos(Rangles[1]), 0, -np.sin(Rangles[1])], [0, 1, 0], [np.sin(Rangles[1]), 0, np.cos(Rangles[1])]])
        r3 = np.array([[np.cos(angles[0]), -np.sin(angles[0]), 0], [np.sin(angles[0]), np.cos(angles[0]), 0], [0, 0, 1]])
        r4 = r3 @ R2 @ R1

        N1 = (A[0, 2] - A[0, 3])
        N2 = (E[0, 1] - C[0, 1])
        f = r4 @ N1
        f = f - N2

        return f[1]

    def M01(self, x0, theta_M, theta_J):
        A, B, C, D, E = self.BaseHexLayer(theta_M, theta_J)

        Rangles = np.zeros(3)

        sol = optimize.root(self.M01_1, x0[0:2], args=(A, C, D, E))
        Rangles[0] = sol.x[0]
        Rangles[1] = sol.x[1]

        sol = optimize.root(self.M01_2, x0[2], args=(Rangles, A, C, E))
        Rangles[2] = sol.x

        R1 = np.array([[1, 0, 0], [0, np.cos(Rangles[0]), np.sin(Rangles[0])], [0, -np.sin(Rangles[0]), np.cos(Rangles[0])]])  # ccw
        R2 = np.array([[np.cos(Rangles[1]), 0, -np.sin(Rangles[1])], [0, 1, 0], [np.sin(Rangles[1]), 0, np.cos(Rangles[1])]])  # ccw
        R3 = np.array([[np.cos(Rangles[2]), -np.sin(Rangles[2]), 0], [np.sin(Rangles[2]), np.cos(Rangles[2]), 0], [0, 0, 1]])  # cw

        matR = R3 @ R2 @ R1

        mat_pack = [A, B, C, D, E]

        return matR, Rangles, mat_pack  # 0th-1st hexagon radial transformation

    def solve_geo_1orbit(self, theta_M: npt.ArrayLike,
                         vtk_out: bool = False, fig_out: bool = True,
                         progbar_out: bool = True):
        progbar_disable = ~progbar_out

        niter = len(theta_M)

        # Import fold angles
        data_th = np.genfromtxt('./lookuptable/theta_MKJS.csv', delimiter=',')
        theta_Mcsv = data_th[:, 0]
        # theta_K = data_th[:, 1]
        theta_J = data_th[:, 2]
        # theta_S = data_th[:, 3]
        # f_MK = interpolate.interp1d(theta_Mcsv, theta_K, kind='cubic')
        f_MJ = interpolate.interp1d(theta_Mcsv, theta_J, kind='cubic')
        # f_MS = interpolate.interp1d(theta_Mcsv, theta_S, kind='cubic')

        theta_J = f_MJ(theta_M)

        A_global = np.zeros((42, 3))
        B_global = np.zeros((42, 3))

        vert_xyz = np.zeros((niter, self.n_node, 3))

        x0 = [0., 0., 0.]

        for i in tqdm(range(niter), disable=progbar_disable):
            matR, Rangles, mat_pack = self.M01(x0, theta_M[i], theta_J[i])
            A, B, C, D, E = mat_pack

            # Update initial guess
            # x0 = Rangles

            for ii in range(0, 6):
                A[1, ii] = E[0, 1] + matR @ (A[0, ii] - A[0, 2])
                B[1, ii] = E[0, 1] + matR @ (B[0, ii] - A[0, 2])
                C[1, ii] = E[0, 1] + matR @ (C[0, ii] - A[0, 2])
                D[1, ii] = E[0, 1] + matR @ (D[0, ii] - A[0, 2])
                E[1, ii] = E[0, 1] + matR @ (E[0, ii] - A[0, 2])
            for ii in range(1, 6):
                for jj in range(0, 5):
                    A[ii + 1, jj + 1] = self.T_6 @ A[ii, jj]
                    A[ii + 1, 0] = self.T_6 @ A[ii, 5]
                    B[ii + 1, jj + 1] = self.T_6 @ B[ii, jj]
                    B[ii + 1, 0] = self.T_6 @ B[ii, 5]
                    C[ii + 1, jj + 1] = self.T_6 @ C[ii, jj]
                    C[ii + 1, 0] = self.T_6 @ C[ii, 5]
                    D[ii + 1, jj + 1] = self.T_6 @ D[ii, jj]
                    D[ii + 1, 0] = self.T_6 @ D[ii, 5]
                    E[ii + 1, jj + 1] = self.T_6 @ E[ii, jj]
                    E[ii + 1, 0] = self.T_6 @ E[ii, 5]

            for ii in range(0, 7):
                for jj in range(0, 6):
                    A_global[jj + 6 * ii] = A[ii, jj]
                    B_global[jj + 6 * ii] = B[ii, jj]

            vert_xyz[i, 0:42, :] = A_global
            vert_xyz[i, 42:48, :] = B_global[0:6]
            vert_xyz[i, 48:54, :] = B_global[[8, 15, 22, 29, 30, 37]]

        if vtk_out:
            self.dir_save = './vtk_resch63_1orbit'
            self.fname_vtk = 'resch63_1orbit'
            if os.path.exists(self.dir_save):
                shutil.rmtree(self.dir_save)
            os.makedirs(self.dir_save)
            print('Create VTK...')
            for i in tqdm(range(niter)):
                self.write_vtk_1orbit(i, vert_xyz[i, :, :])

        if fig_out:
            self.plot_projection(vert_xyz[0], figname='(initial)')
            self.plot_projection(vert_xyz[-1], figname='(final)')

        return vert_xyz

    def write_vtk_1orbit(self, fnum: int, vert_xyz: list):  # , strain: list ):

        with open('%s/%s_%05d.vtk' % (self.dir_save, self.fname_vtk, fnum), 'w') as f:
            f.write('# vtk DataFile Version 3.0\n')
            f.write('Triangulated cylindrical origami (original)\n')
            f.write('ASCII\n')
            num_points = np.size(vert_xyz, axis=0)
            f.write('DATASET POLYDATA\n')
            f.write('POINTS %d float\n' % num_points)
            for i in range(num_points):
                f.write("%f %f %f\n" % (vert_xyz[i, 0], vert_xyz[i, 1], vert_xyz[i, 2]))

            num_dataset = self.n_poly_tri + self.n_poly_hex
            num_datanum = 4 * self.n_poly_tri + 7 * self.n_poly_hex
            f.write('POLYGONS %d %d\n' % (num_dataset, num_datanum))
            # POLYGONS:triangles
            for ip in range(self.n_poly_tri):
                f.write('3 %d %d %d\n' % (self.Polyg_tri[ip, 0], self.Polyg_tri[ip, 1], self.Polyg_tri[ip, 2]))
            # POLYGONS:hexagons
            for ip in range(self.n_poly_hex):
                f.write('6 %d %d %d %d %d %d\n' % (self.Polyg_hex[ip, 0], self.Polyg_hex[ip, 1],
                        self.Polyg_hex[ip, 2], self.Polyg_hex[ip, 3], self.Polyg_hex[ip, 4], self.Polyg_hex[ip, 5]))

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
        for i in range(6):
            ax[0].scatter(vert_xyz[i, 0], vert_xyz[i, 1], s=80, zorder=2,)
        # for i in range(self.n_edge):
        #     ax[0].plot(vert_xyz[self.EdgeConct[i, :], 0],
        #                vert_xyz[self.EdgeConct[i, :], 1],
        #                color='#444444',
        #                linewidth=2,
        #                zorder=1,
        #                )
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
        # for i in range(self.n_edge):
        #     ax[1].plot(vert_xyz[self.EdgeConct[i, :], 0],
        #                vert_xyz[self.EdgeConct[i, :], 2],
        #                color='#444444',
        #                linewidth=2,
        #                zorder=1,
        #                )
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
        # for i in range(self.n_edge):
        #     ax[2].plot(vert_xyz[self.EdgeConct[i, :], 1],
        #                vert_xyz[self.EdgeConct[i, :], 2],
        #                color='#444444',
        #                linewidth=2,
        #                zorder=1,
        #                )
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

    def create_3Dmodel(self, theta_M: npt.ArrayLike,
                       vert_xyz: npt.ArrayLike,
                       dir_save: str = './vtk',
                       fname_vtk: str = 'vtk',):

        # Number of iterations
        niter = len(theta_M)
        # VTK export
        self.dir_save = dir_save
        self.fname_vtk = fname_vtk

        print('Create VTK...')
        # Delete previous data
        if os.path.exists(self.dir_save):
            shutil.rmtree(self.dir_save)
        os.makedirs(self.dir_save)

        # Generate vtk file
        for ii in tqdm(range(niter)):
            # Dummy array for panel color
            strain = (theta_M[ii] - theta_M[0]) / (theta_M[-1] - theta_M[0]) * np.ones(self.n_poly)
            # Call write_vtk to export vtk file at theta_M[ii]
            self.write_vtk(ii, vert_xyz[ii, :, :], strain)

        return

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
            f.write('Triangulated cylindrical origami (original)\n')
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

    def angl2disp(self, theta_M0: float, theta_M: float | npt.ArrayLike,
                  fig_out=False, text_out=False, save_lookuptab: bool = False):
        if isinstance(theta_M, float):
            theta_M = np.array([theta_M, ])
        else:
            theta_M = theta_M

        niter = len(theta_M)
        disp = np.zeros(niter)
        height = np.zeros(niter)

        # Determine coordinates when fully foleded
        vert_xyz = self.solve_geo_1orbit(theta_M=np.array([0, ]),
                                         vtk_out=False, fig_out=False,
                                         progbar_out=False)
        # Extract center hexagon
        vert_xyz_center = vert_xyz[0, :6, :]
        z_center = np.average(vert_xyz_center[:, 2])

        # Determine coordinates at theta_M0
        vert_xyz0 = self.solve_geo_1orbit(theta_M=np.array([theta_M0, ]),
                                          vtk_out=False, fig_out=False,
                                          progbar_out=False)
        # Extract outmost hexagon z coordiante
        z0 = vert_xyz0[0, self.n_node_hex - 6:self.n_node_hex, 2]
        z0 = 0.5 * (max(z0) + min(z0))
        height0 = z_center - z0
        if text_out:
            print('  {0:<32s} : {1:<8f} (deg)'.format('Initial angle', np.degrees(theta_M0)))
            print('  {0:<32s} : {1:<8f}'.format('Initial height', height0))

        vert_xyz = self.solve_geo_1orbit(theta_M=theta_M,
                                         vtk_out=False, fig_out=False,
                                         progbar_out=False)
        zz = vert_xyz[:, self.n_node_hex - 6:self.n_node_hex, 2]
        height = z_center - 0.5 * (np.max(zz, axis=1) + np.min(zz, axis=1))
        disp = height0 - height

        # Save lookup table
        # if save_lookuptab:
        #     print(' Save lookup table')
        #     fpath = './lookuptable/theta_MKJS.csv'
        #     file_exists = os.path.exists(fpath)
        #     if file_exists:
        #         with open(fpath, 'r') as fp:
        #             for count, line in enumerate(fp):
        #                 pass
        #         count = count + 1
        #         file_stat = os.stat(fpath)
        #         if file_stat.st_size >= 1e3 and file_stat.st_size < 1e6:
        #             fsize = file_stat.st_size * 1e-3
        #             funit = 'KB'
        #         elif file_stat.st_size >= 1e6:
        #             fsize = file_stat.st_size * 1e-6
        #             funit = 'MB'
        #         else:
        #             fsize = file_stat.st_size
        #             funit = 'B'
        #         print('  ! File already exists ({0:d} lines; {1:.2f} {2:2s}) !'.format(count, fsize, funit))
        #         value = input('  Overwrite ? [Press Y to overwrite, press elsewhere otherwise]\n')
        #         if value == 'Y' or value == 'y':
        #             print('   Saving lookup table...')
        #             data = np.stack((theta_M, theta_K, theta_J, theta_S), axis=-1)
        #             np.savetxt(fpath, data, delimiter=',')
        #             print('    ...Saved')
        #         else:
        #             print('   File not saved.')
        #             pass
        #     else:
        #         print('   Saving lookup table...')
        #         data = np.stack((theta_M, theta_K, theta_J, theta_S), axis=-1)
        #         np.savetxt(fpath, data, delimiter=',')
        #         print('    ...Saved')

        if fig_out:
            plt.figure('angle-displacement')
            plt.plot(np.degrees(theta_M), disp, label='Displacement $u$')
            plt.plot(np.degrees(theta_M), height, label='Height $h$')
            plt.xlabel('$\\theta_M$')
            plt.ylabel('$u$, $h$')

        return disp, height

    def angl2height(self, theta_M):
        _, height = self.angl2disp(0.0, theta_M, fig_out=False, text_out=False)

        return -height

    def identify_criticalFoldAngle(self,):
        x0 = 0.25 * np.pi
        bnds = ((0, 0.5 * np.pi),)
        options = {'xatol': 1e-8, 'fatol': 1e-8}
        options = {'disp': None,
                   'maxcor': 20,
                   'ftol': 1e-9,
                   'gtol': 1e-8,
                   'eps': 1e-8, }
        soln = optimize.minimize(self.angl2height, x0, bounds=bnds, options=options)
        angl_cr = soln.x[0]

        print(' {0:s} {1:.12f} (deg) = {2:.12f} (rad)'.format('Critical fold angle (thM) is', np.degrees(angl_cr), angl_cr))

        return angl_cr

    def angl2height_res(self, theta_M, h0):
        _, height = self.angl2disp(0.0, theta_M, fig_out=False, text_out=False)

        return height - h0

    def height2angl(self, h0: float, angl_cr: float = 0.74688906):
        x01 = angl_cr * 0.8
        x02 = angl_cr * 1.2
        soln1 = optimize.root(self.angl2height_res, x0=x01, args=(h0))
        soln2 = optimize.root(self.angl2height_res, x0=x02, args=(h0))
        angl_cr1 = soln1.x[0]
        angl_cr2 = soln2.x[0]

        print(' (h0={0:3.2e}) Fold angle (thM) is {1:.12f} (deg) and {2:.12f} (deg)'.format(h0, np.degrees(angl_cr1), np.degrees(angl_cr2)))
        print('                                   {0:.12f} (rad) and {1:.12f} (rad)'.format(angl_cr1, angl_cr2))

        return angl_cr1, angl_cr2

    def analyze_potentialEnergy(self, theta_M0, theta_M,
                                ksp_tor: float = 1,
                                fig_out: bool = True):
        self.ksp_tor = ksp_tor  # torsion spring constant (units:N*m/rad)

        # Import fold angles
        data_th = np.genfromtxt('./lookuptable/theta_MKJS.csv', delimiter=',')
        theta_Mcsv = data_th[:, 0]
        theta_K = data_th[:, 1]
        theta_J = data_th[:, 2]
        theta_S = data_th[:, 3]
        f_MK = interpolate.interp1d(theta_Mcsv, theta_K, kind='cubic')
        f_MJ = interpolate.interp1d(theta_Mcsv, theta_J, kind='cubic')
        f_MS = interpolate.interp1d(theta_Mcsv, theta_S, kind='cubic')

        theta_K0 = f_MK(theta_M0)
        theta_J0 = f_MJ(theta_M0)
        theta_S0 = f_MS(theta_M0)

        theta_K = f_MK(theta_M)
        theta_J = f_MJ(theta_M)
        theta_S = f_MS(theta_M)

        self.n_tor_thM = 24
        self.n_tor_thS = 18
        self.n_tor_thJ = 12
        self.n_tor_thK = 24

        pe_thM = 0.5 * self.Lc * self.n_tor_thM * (2 * theta_M - 2 * theta_M0)**2
        pe_thS = 0.5 * self.Lb * self.n_tor_thS * (2 * theta_S - 2 * theta_S0)**2
        pe_thJ = 0.5 * self.La * self.n_tor_thJ * (theta_J - theta_J0)**2
        pe_thK = 0.5 * self.La * self.n_tor_thK * (2 * theta_K - 2 * theta_K0)**2

        pe_total = self.ksp_tor * (pe_thM + pe_thS + pe_thJ + pe_thK)

        disp, _ = self.angl2disp(theta_M0, theta_M, fig_out=False, text_out=False)

        if fig_out:
            plt.figure('Potential energy (thetaM, normalized)')
            plt.plot(np.degrees(theta_M), pe_total)
            plt.xlabel('$\\theta_M$ (deg)')
            plt.ylabel('Potential energy $U$')

            plt.figure('Potential energy (displacement, normalized)')
            plt.plot(disp * 1e3, pe_total)
            plt.xlabel('Displacement $u$ (mm)')
            plt.ylabel('Potential energy $U$')

        return pe_total

    def analyze_force_displacement(self, theta_M0, theta_M,
                                   ksp_tor: float = 1,
                                   fig_out: bool = True,
                                   save_lookuptab: bool = False,
                                   ):

        self.ksp_tor = ksp_tor  # torsion spring constant (units:N*m/rad)
        pe_total = self.analyze_potentialEnergy(theta_M0, theta_M,
                                                ksp_tor,
                                                fig_out)

        disp, _ = self.angl2disp(theta_M0, theta_M, fig_out=False, text_out=False)

        forc_disp = np.gradient(pe_total, disp, edge_order=2)
        forc_angl = np.gradient(pe_total, theta_M, edge_order=2)
        # pe_thM = np.gradient(pe_total, theta_M, edge_order=2)
        thM_u = np.gradient(theta_M, disp, edge_order=2)
        # forc_disp = pe_thM * thM_u

        forc_disp = forc_disp / (ksp_tor * self.La)

        if save_lookuptab:
            print(' Save lookup table')
            fpath = './lookuptable/forc_disp_angl.csv'
            file_exists = os.path.exists(fpath)
            if file_exists:
                with open(fpath, 'r') as fp:
                    for count, line in enumerate(fp):
                        pass
                count = count + 1
                file_stat = os.stat(fpath)
                if file_stat.st_size >= 1e3 and file_stat.st_size < 1e6:
                    fsize = file_stat.st_size * 1e-3
                    funit = 'KB'
                elif file_stat.st_size >= 1e6:
                    fsize = file_stat.st_size * 1e-6
                    funit = 'MB'
                else:
                    fsize = file_stat.st_size
                    funit = 'B'
                print('  ! File already exists ({0:d} lines; {1:.2f} {2:2s}) !'.format(count, fsize, funit))
                value = input('  Overwrite ? [Press Y to overwrite, press elsewhere otherwise]\n')
                if value == 'Y' or value == 'y':
                    print('   Saving lookup table...')
                    data = np.stack((forc_disp, disp, forc_angl, theta_M), axis=-1)
                    np.savetxt(fpath, data, delimiter=',')
                    print('    ...Saved')
                else:
                    print('   File not saved.')
                    pass
            else:
                print('   Saving lookup table...')
                data = np.stack((forc_disp, disp, forc_angl, theta_M), axis=-1)
                np.savetxt(fpath, data, delimiter=',')
                print('    ...Saved')

        if fig_out:
            plt.figure('Force-angle (normalized)')
            plt.plot(np.degrees(theta_M), forc_angl)
            plt.xlabel('$\\theta_M$ (deg)')
            plt.ylabel('Force $F$')

            plt.figure('Force-displacement (normalized)')
            plt.plot(disp, forc_disp)
            plt.xlabel('Displacement $u$ (m)')
            plt.ylabel('Force $F/(k_0a_0)$')

            plt.figure('Angle-displacement')
            plt.plot(disp * 1e3, np.degrees(theta_M))
            plt.xlabel('Displacement $u$ (mm)')
            plt.ylabel('Angle $\\theta_M$ (deg)')

            plt.figure('Angle-displacement slope')
            plt.plot(disp * 1e3, thM_u)
            plt.xlabel('Displacement $u$ (mm)')
            plt.ylabel('Angle gradient $d\\theta_M/du$ (deg/m)')

        return disp, forc_disp

    def estimate_stiffness(self,):

        return

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
            self.plot_projection(vert_xyz, figname='(init)')

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

    def create_3Dmodel_simplefold(self, theta_M: list | np.ndarray):

        niter = len(theta_M)
        # VTK export
        self.dir_save = './vtk_simplefold'
        self.fname_vtk = 'simplefold'

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

        return


if __name__ == "__main__":

    a = 1.0

    nbin = 100
    theta_M = np.linspace(0., 0.5 * np.pi, nbin)
    # theta_M=np.linspace(np.pi/2, 0, N)
    x0 = np.radians(0.1)

    ROA = ReschOrigamiAnalysis()
    # theta_MKJS = ROA.find_theta_MKJS(theta_M=theta_M,
    #                                  a=a, x0=x0,
    #                                  save_lookuptab=True,
    #                                  )

    theta_M = np.linspace(0.25 * np.pi, np.pi, 500)
    ROA.geo_init_simplefold(La=a, theta_M0=theta_M[0])
    ROA.create_3Dmodel_simplefold(theta_M)

    plt.show()
    exit()
