def main():
    import numpy as np
    import matplotlib.pyplot as plt
    import argparse

    plt.style.use('./common/custom.mplstyle')  # Relative to root directory

    from .lib.miuraori_main import MiuraOriAnalysis
    moa = MiuraOriAnalysis()

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-la', '--lengtha', help='Side length of crease "a" (float, default=1.0)', type=float, default=1.0)
    parser.add_argument('-lb', '--lengthb', help='Side length of crease "b" (float, default=1.0)', type=float, default=1.0)
    parser.add_argument('-ld', '--lengthd', help='Side length of crease "d" (float, default=1.0)', type=float, default=1.0)
    parser.add_argument('-alpha', '--alpha', help='Apex angle "alpha" (float, default=60 deg)', type=float, default=60.0)
    parser.add_argument('-nbin', '--nbin', help='Number of steps', type=int, default=64)
    parser.add_argument('-th0', '--theta0', help='Initial fold angle (float, default=0 deg)', type=float, default=0.0)
    parser.add_argument('-thf', '--thetaf', help='Final fold angle (float, default=90 deg)', type=float, default=90.0)
    parser.add_argument('-zip', '--savezip', help='Option flag to compress vtk files into zip file', action='store_true')
    parser.add_argument('-fig', '--figout', help='Option flag to display plots', action='store_true')

    args = parser.parse_args()
    La = args.lengtha
    Lb = args.lengthb
    Ld = args.lengthb
    alpha = np.radians(args.alpha)
    nbin = args.nbin
    theta0 = np.radians(args.theta0)
    thetaf = np.radians(args.thetaf)
    fig_out = args.figout
    save_zip = args.savezip

    theta_M = np.linspace(theta0, thetaf, nbin)
    moa.geo_init_miura(La=La, Lb=Lb, Ld=Ld,
                       alpha=alpha,
                       theta_M0=theta_M[0],
                       fig_out=fig_out,
                       n_unitx=1, n_unity=1)
    moa.create_3Dmodel(theta_M=theta_M,
                       save_zip=save_zip)

    plt.show()

    return
