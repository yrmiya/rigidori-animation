def main():
    import numpy as np
    import matplotlib.pyplot as plt
    plt.style.use('./common/custom.mplstyle')  # Relative to root directory

    import argparse

    parser = argparse.ArgumentParser(prog='rigidori_animation', description='')
    parser.add_argument('-ori', '--ori_type', help='Type of origami (crease|miura|resch)',
                        choices=['crease', 'miura', 'resch'],
                        type=str, required=True)
    # Universal arguments
    parser.add_argument('-la', '--lengtha', help='Side length of crease "a" (float, default=1.0)', type=float, default=1.0)
    parser.add_argument('-nbin', '--nbin', help='Number of steps', type=int, default=64)
    parser.add_argument('-th0', '--theta0', help='Initial fold angle (float, default=0 deg)', type=float, default=0.0)
    parser.add_argument('-thf', '--thetaf', help='Final fold angle (float, default=90 deg)', type=float, default=90.0)
    parser.add_argument('-zip', '--savezip', help='Option flag to compress vtk files into zip file', action='store_true')
    parser.add_argument('-fig', '--figout', help='Option flag to display plots', action='store_true')
    # Unique to miura ori
    group_miura = parser.add_argument_group('Miura-ori', 'Unique to Miura-ori')
    group_miura.add_argument('-lb', '--lengthb', help='Side length of crease "b" (float, default=1.0)', type=float, default=1.0)
    group_miura.add_argument('-ld', '--lengthd', help='Side length of crease "d" (float, default=1.0)', type=float, default=1.0)
    group_miura.add_argument('-alpha', '--alpha', help='Apex angle "alpha" (float, default=60 deg)', type=float, default=60.0)
    # Unique to resch
    group_resch = parser.add_argument_group('Resch-patterned origami', 'Unique to Resch-patterned origami')
    group_resch.add_argument('-lkup', '--lookuptab',
                             help='Option flag to use lookup table instead of solving fold angle',
                             action='store_false')

    args = parser.parse_args()
    ori_type = args.ori_type
    La = args.lengtha
    nbin = args.nbin
    theta0 = np.radians(args.theta0)
    thetaf = np.radians(args.thetaf)
    fig_out = args.figout
    save_zip = args.savezip

    theta_M = np.linspace(theta0, thetaf, nbin)

    match ori_type:
        case 'crease':
            from .lib.singlecrease import SingleCreaseAnalysis
            sca = SingleCreaseAnalysis()
            sca.geo_init_simplefold(La=La, theta_M0=theta_M[0], fig_out=fig_out)
            sca.create_3Dmodel_simplefold(theta_M, save_zip=save_zip)

        case 'miura':
            Lb = args.lengthb
            Ld = args.lengthb
            alpha = np.radians(args.alpha)
            from .lib.miuraori import MiuraOriAnalysis
            moa = MiuraOriAnalysis()

            moa.geo_init_miura(La=La, Lb=Lb, Ld=Ld,
                               alpha=alpha,
                               theta_M0=theta_M[0],
                               fig_out=fig_out,
                               n_unitx=1, n_unity=1)
            moa.create_3Dmodel(theta_M=theta_M,
                               save_zip=save_zip)
        case 'resch':
            use_lookuptab = args.lookuptab
            from .lib.resch_hex import ReschOrigamiAnalysis
            roa = ReschOrigamiAnalysis()

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
