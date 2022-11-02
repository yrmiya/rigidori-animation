def main():
    import numpy as np
    import matplotlib.pyplot as plt
    import argparse

    plt.style.use('./common/custom.mplstyle')  # Relative to root directory

    from .lib.miuraori_main import MiuraOriAnalysis
    moa = MiuraOriAnalysis()

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-la', '--sidelengtha', help='Side length of crease "a"', type=float, default=1.0)
    parser.add_argument('-lb', '--sidelengthb', help='Side length of crease "b"', type=float, default=1.0)
    parser.add_argument('-nbin', '--nbin', help='Number of steps', type=int, default=64)

    args = parser.parse_args()
    La = args.sidelength
    nbin = args.nbin

    theta_M = np.linspace(0.25 * np.pi, np.pi, nbin)
    moa.geo_init(La=La, theta_M0=theta_M[0])
    moa.create_3Dmodel(theta_M)

    plt.show()

    return
