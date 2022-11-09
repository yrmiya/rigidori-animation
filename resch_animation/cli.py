def main():
    import numpy as np
    import matplotlib.pyplot as plt
    plt.style.use('./common/custom.mplstyle')  # Relative to root directory

    import argparse

    plt.style.use('./common/custom.mplstyle')  # Relative to root directory

    from .lib.resch_hex_static import ReschOrigamiAnalysis
    roa = ReschOrigamiAnalysis()

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-la', '--lengtha', help='Side length of crease "a" (float, default=1.0)', type=float, default=1.0)
    parser.add_argument('-nbin', '--nbin', help='Number of steps', type=int, default=64)
    parser.add_argument('-th0', '--theta0', help='Initial fold angle (float, default=0 deg)', type=float, default=0.0)
    parser.add_argument('-thf', '--thetaf', help='Final fold angle (float, default=90 deg)', type=float, default=90.0)
    parser.add_argument('-zip', '--savezip', help='Option flag to compress vtk files into zip file', action='store_true')
    parser.add_argument('-fig', '--figout', help='Option flag to display plots', action='store_true')
    parser.add_argument('-lkup', '--lookuptab',
                        help='Option flag to use lookup table instead of solving fold angle', action='store_false')

    args = parser.parse_args()
    La = args.lengtha
    nbin = args.nbin
    theta0 = np.radians(args.theta0)
    thetaf = np.radians(args.thetaf)
    fig_out = args.figout
    save_zip = args.savezip
    use_lookuptab = args.lookuptab

    theta_M = np.linspace(theta0, thetaf, nbin)

    if use_lookuptab:
        theta_MKJS = None
    else:
        xx = np.linspace(0.0, 0.5 * np.pi, int(2 * nbin))
        theta_MKJS = roa.find_theta_MKJS(theta_M=xx,
                                         La=La, x0=np.radians(0.1),
                                         save_lookuptab=True,
                                         )

    roa.geo_init(La=La,
                 theta_M0=theta_M[0],
                 theta_MKJS=theta_MKJS,
                 use_lookuptab=use_lookuptab,
                 )

    roa.solve_geo_1orbit(theta_M=theta_M,
                         vtk_out=True,
                         fig_out=fig_out,
                         save_zip=save_zip)

    plt.show()

    return
