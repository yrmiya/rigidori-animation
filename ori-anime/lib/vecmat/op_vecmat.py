import numpy.typing as npt


def calc_facet_area(vert_xyz: npt.ArrayLike,
                    Polyg: npt.ArrayLike,):

    import numpy as np
    n_poly = Polyg.shape[0]
    facet_area = np.zeros(n_poly)
    for ip in range(n_poly):
        vec1 = vert_xyz[Polyg[ip][0], :] - vert_xyz[Polyg[ip][1], :]
        vec2 = vert_xyz[Polyg[ip][2], :] - vert_xyz[Polyg[ip][1], :]
        facet_area[ip] = np.linalg.norm(np.cross(vec1, vec2))

    return facet_area


def calc_rotationMatrix_X(th):
    import numpy as np
    rotT = np.array([[1, 0, 0], [0, np.cos(th), -np.sin(th), 0], [0, np.sin(th), np.cos(th)]])

    return rotT


def calc_rotationMatrix_Y(th):
    import numpy as np
    rotT = np.array([[np.cos(th), 0, -np.sin(th), 0], [0, 1, 0], [np.sin(th), 0, np.cos(th)]])

    return rotT
