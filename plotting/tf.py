import matplotlib.pyplot as plt
import numpy as np


def unwrap(phis, phii, phie):
    phis = (phis - phii + np.pi)%(2*np.pi) - np.pi
    phie = (phie - phii + np.pi)%(2*np.pi) - np.pi
    return phis + phii, phie + phii


def plot_2d_line(rs, zs, re, ze, color=None):
    plt.plot([rs, re], [zs, ze], color=color)


def plot_2d_arc(rs, zs, ri, zi, re, ze, rc, zc, color=None):
    a = np.linalg.norm([ri - rc, zi - zc])
    assert(np.isclose(np.linalg.norm([rs - rc, zs - zc]), a) and np.isclose(np.linalg.norm([re - rc, ze - zc]), a))
    ths = np.arctan2(zs - zc, rs - rc)
    thi = np.arctan2(zi - zc, ri - rc)
    the = np.arctan2(ze - zc, re - rc)
    ths, the = unwrap(ths, thi, the)

    n = 33
    th = np.linspace(ths, the, n)
    plt.plot(rc + a*np.cos(th), zc + a*np.sin(th), color=color)


def plot_2d_conductor(conductor, color=None):
    elements = conductor.elements
    start_points = elements.start_points
    intermediate_points = elements.intermediate_points
    end_points = elements.end_points
    centres = elements.centres
    for typ, rs, zs, ri, zi, re, ze, rc, zc in zip(elements.types, start_points.r, start_points.z, intermediate_points.r, intermediate_points.z, end_points.r, end_points.z, centres.r, centres.z):
        if typ == 1:
            plot_2d_line(rs, zs, re, ze, color=color)
        elif typ == 2:
            plot_2d_arc(rs, zs, ri, zi, re, ze, rc, zc, color=color)
        else:
            raise RuntimeError("Not implemented")


def plot_2d_coil(coil, color=None):
    for conductor in coil.conductor:
        plot_2d_conductor(conductor, color=color)


def plot_2d(tf):
    for coil in tf.coil:
        plot_2d_coil(coil, color="silver")


def plot_tf(tf):
    n = len(tf.coil)
    r = int(np.sqrt(n))
    c = (n + r - 1)//r

    fig, axs = plt.subplots(r, c, sharex=True)
    for ax, coil in zip(np.atleast_1d(axs).flatten(), tf.coil):
        ax.plot(tf.time, coil.current.data)
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
        plot_tf(db.get("tf"))

    plt.show()
