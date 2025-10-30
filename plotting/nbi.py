import matplotlib.pyplot as plt
import numpy as np


def unit_vector(unit_vector):
    return np.array([unit_vector.x, unit_vector.y, unit_vector.z])


def rphiz2xyz(rphiz):
    return np.array([rphiz.r*np.cos(rphiz.phi), rphiz.r*np.sin(rphiz.phi), rphiz.z])


def xyz2rz(xyz):
    return np.array([np.linalg.norm(xyz[:2], axis=0), xyz[2]])


def plot_2d_rectangle(r, x1, x2, width, height, color=None):
    plt.plot(*xyz2rz(np.array([r + 0.5*x1*width + 0.5*x2*height, r - 0.5*x1*width + 0.5*x2*height]).T), color=color)
    plt.plot(*xyz2rz(np.array([r - 0.5*x1*width + 0.5*x2*height, r - 0.5*x1*width - 0.5*x2*height]).T), color=color)
    plt.plot(*xyz2rz(np.array([r - 0.5*x1*width - 0.5*x2*height, r + 0.5*x1*width - 0.5*x2*height]).T), color=color)
    plt.plot(*xyz2rz(np.array([r + 0.5*x1*width - 0.5*x2*height, r + 0.5*x1*width + 0.5*x2*height]).T), color=color)


def plot_2d_unit(unit, color=None):
    aperture = unit.aperture[0]
    beamlets_group = unit.beamlets_group[0]

    x1_unit_vector = np.array([aperture.x1_unit_vector.x, aperture.x1_unit_vector.y, aperture.x1_unit_vector.z])
    x2_unit_vector = np.array([aperture.x2_unit_vector.x, aperture.x2_unit_vector.y, aperture.x2_unit_vector.z])
    x3_unit_vector = np.array([aperture.x3_unit_vector.x, aperture.x3_unit_vector.y, aperture.x3_unit_vector.z])

    position = beamlets_group.position
    xsc = rphiz2xyz(position)
    xap = rphiz2xyz(aperture.centre)

    lbsctan = (position.r**2 - beamlets_group.tangency_radius**2)**0.5/np.cos(beamlets_group.angle)
    tangency_point = xsc + x3_unit_vector*lbsctan

    plot_2d_rectangle(
        tangency_point,
        x1_unit_vector,
        x2_unit_vector,
        beamlets_group.width_horizontal*(1 - lbsctan/beamlets_group.focus.focal_length_horizontal),
        beamlets_group.width_vertical*(1 - lbsctan/beamlets_group.focus.focal_length_vertical),
        color,
    )


def plot_2d(nbi):
    for unit, color in zip(nbi.unit, ['r', 'g', 'b', 'm', 'k', 'c']):
        plot_2d_unit(unit, color)


def plot_3d_rectangle(r, x1, x2, width, height, color=None):
    plt.plot(*np.array([r + 0.5*x1*width + 0.5*x2*height, r - 0.5*x1*width + 0.5*x2*height]).T, color=color)
    plt.plot(*np.array([r - 0.5*x1*width + 0.5*x2*height, r - 0.5*x1*width - 0.5*x2*height]).T, color=color)
    plt.plot(*np.array([r - 0.5*x1*width - 0.5*x2*height, r + 0.5*x1*width - 0.5*x2*height]).T, color=color)
    plt.plot(*np.array([r + 0.5*x1*width - 0.5*x2*height, r + 0.5*x1*width + 0.5*x2*height]).T, color=color)


def plot_3d_unit(unit, color=None):
    # Assume one beamlets_group, one aperture
    # Assume x1, x2, x3 unit vectors are shared for aperture and beamlets_group
    aperture = unit.aperture[0]
    beamlets_group = unit.beamlets_group[0]

    x1_unit_vector = np.array([aperture.x1_unit_vector.x, aperture.x1_unit_vector.y, aperture.x1_unit_vector.z])
    x2_unit_vector = np.array([aperture.x2_unit_vector.x, aperture.x2_unit_vector.y, aperture.x2_unit_vector.z])
    x3_unit_vector = np.array([aperture.x3_unit_vector.x, aperture.x3_unit_vector.y, aperture.x3_unit_vector.z])

    position = beamlets_group.position
    xsc = rphiz2xyz(position)
    xap = rphiz2xyz(aperture.centre)

    lbsctan = (position.r**2 - beamlets_group.tangency_radius**2)**0.5/np.cos(beamlets_group.angle)
    lbscap = (xap - xsc).dot(x3_unit_vector)
    tangency_point = xsc + x3_unit_vector*lbsctan
    plt.plot(*np.vstack([xsc, tangency_point]).T, f"{color}:")

    plot_3d_rectangle(
        xsc,
        x1_unit_vector,
        x2_unit_vector,
        beamlets_group.width_horizontal,
        beamlets_group.width_vertical,
        color,
    )

    plot_3d_rectangle(
        xap,
        x1_unit_vector,
        x2_unit_vector,
        aperture.x1_width,
        aperture.x2_width,
        color,
    )

    plot_3d_rectangle(
        tangency_point,
        x1_unit_vector,
        x2_unit_vector,
        beamlets_group.width_horizontal*(1 - lbsctan/beamlets_group.focus.focal_length_horizontal),
        beamlets_group.width_vertical*(1 - lbsctan/beamlets_group.focus.focal_length_vertical),
        color,
    )

    # Check these values with NUBEAM parameters
    """
    srtcen = beamlets_group.direction*beamlets_group.tangency_radius
    lbsctan = (position.r**2 - beamlets_group.tangency_radius**2)**0.5/np.cos(beamlets_group.angle)
    zbsc = position.z
    phibsc = rad2deg*position.phi
    lbscap = (xap - xsc).dot(x3_unit_vector)
    zbap = zbsc + x3_unit_vector[2]*lbscap
    b_halfwidth = 0.5*beamlets_group.width_horizontal
    b_halfheight = 0.5*beamlets_group.width_vertical
    b_hfocal_length = beamlets_group.focus.focal_length_horizontal
    b_vfocal_length = beamlets_group.focus.focal_length_vertical
    b_hdivergence = rad2deg*beamlets_group.divergence_component[0].horizontal
    b_vdivergence = rad2deg*beamlets_group.divergence_component[0].vertical
    ap_halfwidth = 0.5*aperture.x1_width
    ap_halfheight = 0.5*aperture.x2_width
    ap_horiz_offset = (xap - xsc).dot(x1_unit_vector)
    ap_vert_offset = -(xap - xsc).dot(x2_unit_vector)
    """


def plot_3d(nbi):
    for unit, color in zip(nbi.unit, ['r', 'g', 'b', 'm', 'k', 'c']):
        plot_3d_unit(unit, color)


def plot_nbi(nbi):
    n = len(nbi.unit)
    r = int(np.sqrt(n))
    c = (n + r - 1)//r

    fig, axs = plt.subplots(r, c, sharex=True)
    for ax, unit in zip(np.atleast_1d(axs).flatten(), nbi.unit):
        ax.set_ylabel(r"Power [$W$], Energy [$0.1 V$]")
        ax.plot(nbi.time, unit.power_launched.data, color='r')
        ax.plot(nbi.time, 10.0*unit.energy.data, color='b')

        tx = ax.twinx()
        tx.set_ylabel(r"Beam current fraction")
        tx.plot(nbi.time, unit.beam_current_fraction.data.T)
        tx.set_ylim([0.0, 1.0])

        ax.set_title(unit.identifier)


if __name__ == "__main__":
    import imas
    import sys

    try:
        path = sys.argv[1]
    except IndexError:
        print(f"Usage: {sys.argv[0]} $PATH")
        exit()

    with imas.DBEntry(f"imas:hdf5?path={path}", 'r') as db:
        plot_nbi(db.get("nbi"))

    plt.show()
