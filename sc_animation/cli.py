def main():
    import numpy as np
    import matplotlib.pyplot as plt
    import argparse

    plt.style.use('./common/custom.mplstyle')  # Relative to root directory

    from .lib.singlecrease_main import SingleCreaseAnalysis
    sca = SingleCreaseAnalysis()

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-sl', '--sidelength', help='Side length (float)', type=float, default=1.0)
    parser.add_argument('-nbin', '--nbin', help='Number of steps (int)', type=int, default=64)
    parser.add_argument('-th0', '--theta0', help='Initial fold angle (float)', type=float, default=0.25 * np.pi)
    parser.add_argument('-thf', '--thetaf', help='Final fold angle (float)', type=float, default=np.pi)
    parser.add_argument('-zip', '--savezip', help='Option to compress vtk files into zip file (True|False)', type=bool, default=False)
    parser.add_argument('-fig', '--figout', help='Option to display plots (True|False)', type=bool, default=False)

    args = parser.parse_args()
    La = args.sidelength
    nbin = args.nbin
    theta0 = args.theta0
    thetaf = args.thetaf
    fig_out = args.figout
    save_zip = args.savezip

    theta_M = np.linspace(theta0, thetaf, nbin)
    sca.geo_init_simplefold(La=La, theta_M0=theta_M[0], fig_out=fig_out)
    sca.create_3Dmodel_simplefold(theta_M, save_zip=save_zip)

    plt.show()

    return
