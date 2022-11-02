def main():
    import numpy as np
    import matplotlib.pyplot as plt
    plt.style.use('./common/custom.mplstyle')  # Relative to root directory

    from .lib.resch_hex_static import ReschOrigamiAnalysis
    from .lib.cyclic_main import CyclicTestAnalysis
    from .lib.read_experiment import QuasiStaticExperimentAnalysis
    ROA = ReschOrigamiAnalysis()
    CTA = CyclicTestAnalysis()
    QSEA = QuasiStaticExperimentAnalysis()

    La = 1.0
    La = 30.0e-3
    La = 30.

    # =========================================================
    #   Find theta_K, theta_J, and theta_S
    #   as a function of theta_M
    #   comment out to skip
    # =========================================================
    # nbin = 256
    # theta_M = np.linspace(0., 0.5 * np.pi, nbin)
    # x0 = np.radians(0.1)

    # theta_MKJS = ROA.find_theta_MKJS(theta_M=theta_M,
    #                                  Lc=Lc, x0=x0,
    #                                  save_lookuptab=True,
    #                                  )

    # =========================================================
    #   Initialize geometry at given theta_M0
    #   comment out to skip
    # =========================================================
    theta_M0 = 0.25 * np.pi
    ROA.geo_init(La=La, n_orbit=1, theta_M0=theta_M0)
    # exit()

    # =========================================================
    #   Solve folding behavior for given series of theta_M
    #   comment out to skip
    # =========================================================
    nbin = 256
    theta_M = np.linspace(0., np.pi / 12., nbin)
    # ROA.solve_geo_1orbit(theta_M=theta_M, vtk_out=False, fig_out=True)
    # plt.show()
    # exit()

    # =========================================================
    #   Identify critical fold angle theta_M
    #   comment out to skip
    # =========================================================
    # thM_cr = ROA.identify_criticalFoldAngle()
    # Critical fold angle (thM) is 42.79359091 (deg) = 0.74688906 (rad)

    # =========================================================
    #   Fold angle to displcaement + height conversion
    #   comment out to skip
    # =========================================================
    nbin = 512
    thM_cr = 0.746889060157
    theta_M0 = np.pi / 2.
    theta_M1 = 0.0
    theta_M = np.linspace(theta_M0, theta_M1, nbin)
    # ROA.angl2disp(theta_M0, theta_M, fig_out=True)

    # plt.show()
    # exit()

    # =========================================================
    #   Identify fold angle theta_M at given initial height h0
    #   comment out to skip
    # =========================================================
    h0 = 10
    h0 = 12
    # thM_soln = ROA.height2angl(h0)
    h0 = 15
    # thM_soln = ROA.height2angl(h0)

    # plt.show()
    # exit()

    # =========================================================
    #   Energy-displacemenet
    #   comment out to skip
    # =========================================================
    nbin = 128
    theta_M0 = np.pi / 2.
    theta_M1 = 0.0
    theta_M0 = 0.0
    theta_M1 = np.pi / 2.
    theta_M = np.linspace(theta_M0, theta_M1, nbin)
    # ROA.analyze_potentialEnergy(theta_M0, theta_M, )

    # =========================================================
    #   Force-displacemenet
    #   comment out to skip
    # =========================================================
    nbin = 128
    nbin = 1024
    thM_cr = 0.746889060157  # (rad)
    ksp_tor = 1.0

    # theta_M0 = thM_soln[0]
    theta_M0 = 0.6308865206991643  # h0=15
    if theta_M0 > thM_cr:
        theta_M1 = 0.5 * np.pi
    elif theta_M0 < thM_cr:
        theta_M1 = 0.0
    theta_M = np.linspace(theta_M0, theta_M1, nbin)
    ROA.analyze_force_displacement(theta_M0, theta_M, fig_out=True)

    theta_M0 = np.radians(23.485184170197)
    if theta_M0 > thM_cr:
        theta_M1 = 0.5 * np.pi
    elif theta_M0 < thM_cr:
        theta_M1 = 0.0
    theta_M = np.linspace(theta_M0, theta_M1, nbin)
    # ROA.analyze_force_displacement(theta_M0, theta_M, fig_out=True)

    theta_M0 = np.radians(62.548873947182)
    if theta_M0 > thM_cr:
        theta_M1 = 0.5 * np.pi
    elif theta_M0 < thM_cr:
        theta_M1 = 0.0
    theta_M = np.linspace(theta_M0, theta_M1, nbin)
    # ROA.analyze_force_displacement(theta_M0, theta_M, fig_out=True)

    # =========================================================
    #   Analyze experimental results
    #   comment out to skip
    # =========================================================
    # ********************************
    #   Single compression analysis
    # ********************************
    val = True
    val = False
    if val:
        filename1 = './data/experiment/20220914_Resch63_O1C30_P'
        filename2 = '_N10'
        id_file = [[1, ], ]
        id_prototype = [1, ]
        filename1 = './data/experiment/20220914_Resch63_O1C30_P'
        filename2 = ''
        id_file = [[1, 2, 3, ], ]
        id_prototype = [1, ]

        disp1, f_ave1, f_err1, g_ave1 = QSEA.plot_experiment(ftype='cycl',  # ftype: 'none' or 'cycl'
                                                             filename1=filename1,
                                                             filename2=filename2,
                                                             id_file=id_file,
                                                             id_prototype=id_prototype,
                                                             ex_type='comp',
                                                             color='#f08687',
                                                             fig_out=True,
                                                             data_trim=[None, None],
                                                             )

    # *****************************
    #   Cyclic loading analysis
    # *****************************
    val = True
    val = False
    if val:
        # index_fd_plot: Cycle(s) to be plotted on force-displacement curve
        fpath = '../../data/experiment/kresling_2017intern'
        date = '20170714'
        filename = 'TCO_H35C36Q70_A'
        figname = 'monostable(A,N200)'
        sample_i = 1
        sample_f = 5
        index_fd_plot = [1, 2, 5, 10, 50, 100]

        fpath = './data/experiment'
        date = '20221005'
        filename = 'Resch63_O1C30N50_PX'
        figname = 'Resch63'
        sample_i = 3
        sample_f = 5
        index_fd_plot = [1, 2, 5, 10, 50]

        CTA.PlotData(fpath=fpath,
                     date=date,
                     filename=filename,
                     sample_i=sample_i,
                     sample_f=sample_f,
                     # sample_f=23,
                     figname=figname,
                     color='#f08687',
                     bbox_to_anchor=(1.05, 0.4),
                     stability='mono',
                     full_fd=True,
                     index_fd_plot=index_fd_plot
                     )

    # =========================================================
    #   Example code to create folding animation of
    #   single crease
    #   comment out to skip
    # =========================================================
    # theta_M = np.linspace(0.25 * np.pi, np.pi, 500)
    # ROA.geo_init_simplefold(Lc=a, theta_M0=theta_M[0])
    # ROA.create_3Dmodel_simplefold(theta_M)

    plt.show()

    return
