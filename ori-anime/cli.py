def main():
    import argparse

    parser = argparse.ArgumentParser(prog='rigidori_animation', description='')

    subparsers = parser.add_subparsers(help='Sub-commands')
    parser_main = subparsers.add_parser('run', help='Execute tool')
    parser_main.set_defaults(handler=default_run)
    parser_clean = subparsers.add_parser('clean', help='Remove previous results')
    parser_clean.set_defaults(handler=clean_all)

    parser_main.add_argument('-ori', '--ori_type', help='Type of origami (crease|miura|resch)',
                             choices=['crease', 'miura', 'resch', 'startuck', 'kres'],
                             type=str, required=True)
    # Universal arguments
    parser_main.add_argument('-la', '--lengtha', help='Side length of crease "a" (float, default=1.0)', type=float, default=1.0)
    parser_main.add_argument('-nbin', '--nbin', help='Number of steps', type=int, default=64)
    parser_main.add_argument('-th0', '--theta0', help='Initial fold angle (float, default=0 deg)', type=float, default=0.0)
    parser_main.add_argument('-thf', '--thetaf', help='Final fold angle (float, default=90 deg)', type=float, default=90.0)
    parser_main.add_argument('-zip', '--savezip', help='Option flag to compress vtk files into zip file', action='store_true')
    parser_main.add_argument('-fig', '--figout', help='Option flag to display plots', action='store_true')
    parser_main.add_argument('-clean', '--clean', help='Option flag to remove vtk files after creating zip file', action='store_true')
    # Unique to miura ori
    group_miura = parser_main.add_argument_group('Miura-ori', 'Unique to Miura-ori')
    group_miura.add_argument('-lb', '--lengthb', help='Side length of crease "b" (float, default=1.0)', type=float, default=1.0)
    group_miura.add_argument('-ld', '--lengthd', help='Side length of crease "d" (float, default=1.0)', type=float, default=1.0)
    group_miura.add_argument('-alpha', '--alpha', help='Apex angle "alpha" (float, default=60 deg)', type=float, default=60.0)
    group_miura.add_argument('-unitx', '--unitx', help='Number of unit in x direction', type=int, default=1)
    # Unique to Kresling
    group_kresling = parser_main.add_argument_group('Kresling', 'Unique to Kresling-ori')
    group_kresling.add_argument('-h0', '--h0', help='Initial height (float, default=1)', type=float, default=1.0)
    group_kresling.add_argument('-ph0', '--ph0', help='Initial rotation (float, default=60 deg)', type=float, default=60.0)
    group_kresling.add_argument('-unit', '--unit', help='Number of unit in x direction', type=int, default=1)
    group_kresling.add_argument('-side', '--side', help='Number of side of polygonal section', type=int, default=6)
    # Unique to resch
    group_resch = parser_main.add_argument_group('Resch-patterned origami', 'Unique to Resch-patterned origami')
    group_resch.add_argument('-lkup', '--lookuptab',
                             help='Option flag to use lookup table instead of solving fold angle',
                             action='store_false')
    # Unique to startuck
    group_star = parser_main.add_argument_group('Startuck fold', 'Unique to Startuck fold')
    group_star.add_argument('-st_type', '--startuck_type',
                            help='Startuck profile',
                            choices=['tri'],
                            type=str,
                            default='tri')
    group_star.add_argument('-bi', '--bistartuck',
                            help='Export two startuck tessellation',
                            action='store_true')

    args = parser.parse_args()
    if hasattr(args, "handler"):
        args.handler(args)
    else:
        parser.print_help()

    return


def clean_all(args):
    import glob
    import os
    import shutil
    fpath1 = './vtk_*.zip'
    dfile1 = glob.glob(fpath1)
    if len(dfile1) > 0:
        for i in range(len(dfile1)):
            os.remove(path=dfile1[i])
        fpath2 = './vtk_*'
        dfile2 = glob.glob(fpath2)
        if len(dfile2) > 0:
            for i in range(len(dfile2)):
                shutil.rmtree(path=dfile2[i])
        print('Cleaned all previous results')
    else:
        fpath2 = './vtk_*'
        dfile2 = glob.glob(fpath2)
        if len(dfile2) > 0:
            for i in range(len(dfile2)):
                shutil.rmtree(path=dfile2[i])
            print('Cleaned all previous results')
        else:
            print('No previous results to remove')
    return


def default_run(args):
    import numpy as np
    import matplotlib.pyplot as plt
    plt.style.use('common/custom.mplstyle')  # Relative to root directory

    ori_type = args.ori_type
    La = args.lengtha
    nbin = args.nbin
    theta0 = np.radians(args.theta0)
    thetaf = np.radians(args.thetaf)
    fig_out = args.figout
    save_zip = args.savezip
    clean_vtk = args.clean

    theta_M = np.linspace(theta0, thetaf, nbin)

    match ori_type:
        case 'crease':
            dir_save = './vtk_singlecrease'
            fname_vtk = 'singlecrease'

            from .lib.singlecrease import SingleCreaseAnalysis
            sca = SingleCreaseAnalysis(dir_save, fname_vtk)
            sca.geo_init_simplefold(La=La, theta_M0=theta_M[0], fig_out=fig_out)
            sca.create_3Dmodel_simplefold(theta_M, save_zip=save_zip)

        case 'miura':
            Lb = args.lengthb
            Ld = args.lengthb
            alpha = np.radians(args.alpha)
            n_unitx = args.unitx

            dir_save = './vtk_miura'
            fname_vtk = 'miura'

            from .lib.miuraori import MiuraOriAnalysis
            moa = MiuraOriAnalysis(dir_save, fname_vtk)

            moa.geo_init_miura(La=La, Lb=Lb, Ld=Ld,
                               alpha=alpha,
                               theta_M0=theta_M[0],
                               fig_out=fig_out,
                               n_unitx=n_unitx, n_unity=1)
            moa.create_3Dmodel(theta_M=theta_M,
                               save_zip=save_zip)

        case 'kres':
            n_unit = args.unit
            n_side = args.side
            Lc = La
            h0 = args.h0
            ph0 = np.radians(args.ph0)

            dir_save = './vtk_kresling'
            fname_vtk = 'kresling'

            from .lib.kresling import KreslingAnalysis
            kres = KreslingAnalysis(dir_save, fname_vtk)

            kres.geo_init_kres(n_unit=n_unit,
                               n_side=n_side,
                               Lc=Lc,
                               h0=h0,
                               ph0=ph0)
            disp = np.linspace(0, h0, nbin)

        case 'resch':
            use_lookuptab = args.lookuptab

            dir_save = './vtk_resch63_1orbit'
            fname_vtk = 'resch63_1orbit'

            from .lib.resch_hex import ReschOrigamiAnalysis
            roa = ReschOrigamiAnalysis(dir_save, fname_vtk)

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

        case 'startuck':
            st_type = args.startuck_type
            from .lib.startuck import StartuckAnalysisTool
            star = StartuckAnalysisTool()

            star.geo_init_startuck(La=La,
                                   theta_M0=theta_M[0],
                                   st_type=st_type,
                                   fig_out=fig_out)

            star.create_3Dmodel(theta_M=theta_M, save_zip=save_zip)
            bistartuck = args.bistartuck
            if bistartuck:
                star.create_3Dmodel_unit(theta_M=theta_M, save_zip=save_zip)

    if save_zip and clean_vtk:
        import shutil
        shutil.rmtree(path=dir_save)

    plt.show()
    return
