def main():
    import numpy as np
    import matplotlib.pyplot as plt
    import argparse

    plt.style.use('./common/custom.mplstyle')  # Relative to root directory

    from .lib.singlecrease_main import SingleCreaseAnalysis
    sca = SingleCreaseAnalysis()

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-sa', '--sidelength', help='Side length', type=float, default=1.0)
    parser.add_argument('-nbin', '--nbin', help='Number of steps', type=int, default=64)

    args = parser.parse_args()
    La = args.sidelength
    nbin = args.nbin

    theta_M = np.linspace(0.25 * np.pi, np.pi, nbin)
    sca.geo_init_simplefold(La=La, theta_M0=theta_M[0])
    sca.create_3Dmodel_simplefold(theta_M)

    plt.show()

    return
