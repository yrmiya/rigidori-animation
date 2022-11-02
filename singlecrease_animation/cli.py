def main():
    import numpy as np
    import matplotlib.pyplot as plt
    plt.style.use('./common/custom.mplstyle')  # Relative to root directory

    from .lib.singlecrease_main import SingleCreaseAnalysis
    sca = SingleCreaseAnalysis()

    La = 1.0

    # =========================================================
    #   Example code to create folding animation of
    #   single crease
    #   comment out to skip
    # =========================================================
    nbin = 512
    theta_M = np.linspace(0.25 * np.pi, np.pi, nbin)
    sca.geo_init_simplefold(La=La, theta_M0=theta_M[0])
    sca.create_3Dmodel_simplefold(theta_M)

    plt.show()

    return
