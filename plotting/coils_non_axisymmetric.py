import matplotlib.pyplot as plt
import numpy as np
import itertools


colors = itertools.cycle(['r', 'g', 'b', 'm', 'k', 'c'])


def rphiz2xyz(r, phi, z):
    return np.array([r*np.cos(phi), r*np.sin(phi), z])


def unwrap(phis, phii, phie):
    phis = (phis - phii + np.pi)%(2*np.pi) - np.pi
    phie = (phie - phii + np.pi)%(2*np.pi) - np.pi
    return phis + phii, phie + phii


def plot_3d_line(rs, phis, zs, re, phie, ze, color=None):
    plt.plot(*np.vstack([rphiz2xyz(rs, phis, zs), rphiz2xyz(re, phie, ze)]).T, color=color)


def plot_3d_arc(rs, phis, zs, ri, phii, zi, re, phie, ze, color=None):
    assert(rs == re)
    assert(zs == ze)

    n = 33
    phis, phie = unwrap(phis, phii, phie)
    phi = np.linspace(phis, phie, n)
    for phis, phie in zip(phi[:-1], phi[1:]):
        plt.plot(*np.vstack([rphiz2xyz(rs, phis, zs), rphiz2xyz(re, phie, ze)]).T, color=color)


def plot_3d_conductor(conductor, color=None):
    s = conductor.elements.start_points
    i = conductor.elements.intermediate_points
    e = conductor.elements.end_points
    #c = conductor.elements.centres
    for typ, rs, phis, zs, ri, phii, zi, re, phie, ze in zip(conductor.elements.types, s.r, s.phi, s.z, i.r, i.phi, i.z, e.r, e.phi, e.z):
        if typ == 1:
            plot_3d_line(rs, phis, zs, re, phie, ze, color)
        elif typ == 2:
            # Assume rc = 0
            plot_3d_arc(rs, phis, zs, ri, phii, zi, re, phie, ze, color)
        else:
            raise RuntimeError("Not implemented")


def plot_3d_coil(coil, color):
    for conductor in coil.conductor:
        plot_3d_conductor(conductor, color)


def plot_3d(coils_non_axisymmetric):
    for coil, color in zip(coils_non_axisymmetric.coil, colors):
        plot_3d_coil(coil, color)


def plot_coils_non_axisymmetric(coils_non_axisymmetric):
    n = len(coils_non_axisymmetric.coil)
    r = int(np.sqrt(n))
    c = (n + r - 1)//r

    fig, axs = plt.subplots(r, c, sharex=True)
    for ax, coil in zip(np.atleast_1d(axs).flatten(), coils_non_axisymmetric.coil):
        ax.plot(coils_non_axisymmetric.time, coil.current.data)
        ax.set_title(coil.identifier)


if __name__ == "__main__":
    import imas
    import sys

    try:
        path = sys.argv[1]
    except IndexError:
        print(f"Usage: {sys.argv[0]} $PATH")
        exit()

    with imas.DBEntry(f"imas:hdf5?path={path}", 'r') as db:
        plot_coils_non_axisymmetric(db.get("coils_non_axisymmetric"))

    plt.show()
