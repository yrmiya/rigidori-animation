import numpy.typing as npt


def plot_projection(vert_xyz: npt.ArrayLike,
                    Polyg: npt.ArrayLike,
                    EdgeConct: npt.ArrayLike,
                    figname: str = ''):
    '''
    Plot xy-plane, xz-plane, and yz-plane projection of the origami

    Args:
        vert_xyz (npt.ArrayLike): Nodal coordinates
        Polyg (npt.ArrayLike): Polygon-node connectivity
        EdgeConct (npt.ArrayLike): Element-node connectivity
        figname (str, optional): Figure name. Defaults to ''.
    '''
    
    import numpy as np
    import matplotlib.pyplot as plt

    n_node = vert_xyz.shape[0]
    n_poly = Polyg.shape[0]
    n_edge = EdgeConct.shape[0]

    fig, ax = plt.subplots(1, 3, num='xyz projection %s' % (figname), figsize=(18, 6))
    # xy plane
    for i in range(n_node):
        ax[0].scatter(vert_xyz[i, 0], vert_xyz[i, 1], s=80, zorder=2,)
    for i in range(n_edge):
        ax[0].plot(vert_xyz[EdgeConct[i, :], 0],
                   vert_xyz[EdgeConct[i, :], 1],
                   color='#444444',
                   linewidth=2,
                   zorder=1,
                   )
    for i in range(n_poly):
        ax[0].fill(vert_xyz[Polyg[i][:], 0],
                   vert_xyz[Polyg[i][:], 1],
                   alpha=0.2,
                   linewidth=2,
                   zorder=0,
                   )
    ax[0].set_xlabel('$x$')
    ax[0].set_ylabel('$y$')
    ax[0].set_aspect('equal', 'box')
    # xz plane
    for i in range(n_node):
        ax[1].scatter(vert_xyz[i, 0], vert_xyz[i, 2], s=80, zorder=2,)
    for i in range(n_edge):
        ax[1].plot(vert_xyz[EdgeConct[i, :], 0],
                   vert_xyz[EdgeConct[i, :], 2],
                   color='#444444',
                   linewidth=2,
                   zorder=1,
                   )
    for i in range(n_poly):
        ax[1].fill(vert_xyz[Polyg[i][:], 0],
                   vert_xyz[Polyg[i][:], 2],
                   alpha=0.2,
                   linewidth=2,
                   zorder=0,
                   )
    ax[1].set_xlabel('$x$')
    ax[1].set_ylabel('$z$')
    ax[1].set_aspect('equal', 'box')
    # yz plane
    for i in range(n_node):
        ax[2].scatter(vert_xyz[i, 1], vert_xyz[i, 2], s=80, zorder=2,)
    for i in range(n_edge):
        ax[2].plot(vert_xyz[EdgeConct[i, :], 1],
                   vert_xyz[EdgeConct[i, :], 2],
                   color='#444444',
                   linewidth=2,
                   zorder=1,
                   )
    for i in range(n_poly):
        ax[2].fill(vert_xyz[Polyg[i][:], 1],
                   vert_xyz[Polyg[i][:], 2],
                   alpha=0.2,
                   linewidth=2,
                   zorder=0,
                   )
    ax[2].set_xlabel('$y$')
    ax[2].set_ylabel('$z$')
    ax[2].set_aspect('equal', 'box')
    fig.tight_layout()

    return
