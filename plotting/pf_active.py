import matplotlib.pyplot as plt
import numpy as np


def plot_2d_coil(coil):
    rl = min(element.geometry.rectangle.r - 0.5*element.geometry.rectangle.width for element in coil.element)
    rr = max(element.geometry.rectangle.r + 0.5*element.geometry.rectangle.width for element in coil.element)
    zl = min(element.geometry.rectangle.z - 0.5*element.geometry.rectangle.height for element in coil.element)
    zr = max(element.geometry.rectangle.z + 0.5*element.geometry.rectangle.height for element in coil.element)

    r = np.array([rr, rl, rl, rr, rr])
    z = np.array([zr, zr, zl, zl, zr])
    plt.plot(r, z, color='k')


def plot_2d(pf_active):
    for coil in pf_active.coil:
        plot_2d_coil(coil)


def plot_3d_coil(coil, nphi=65):
    phi = np.linspace(0, 2*np.pi, nphi)

    rl = min(element.geometry.rectangle.r - 0.5*element.geometry.rectangle.width for element in coil.element)
    rr = max(element.geometry.rectangle.r + 0.5*element.geometry.rectangle.width for element in coil.element)
    zl = min(element.geometry.rectangle.z - 0.5*element.geometry.rectangle.height for element in coil.element)
    zr = max(element.geometry.rectangle.z + 0.5*element.geometry.rectangle.height for element in coil.element)

    r = np.array([rr, rl, rl, rr, rr])
    z = np.array([zr, zr, zl, zl, zr])
    plt.gca().plot_surface(
        np.outer(r, np.cos(phi)),
        np.outer(r, np.sin(phi)),
        np.array([z]*nphi).T,
        color='k',
        alpha=0.5,
    )


def plot_3d(pf_active, nphi=65):
    for coil in pf_active.coil:
        plot_3d_coil(coil, nphi)


def plot_pf_active(pf_active):
    n = len(pf_active.coil)
    r = int(np.sqrt(n))
    c = (n + r - 1)//r

    fig, axs = plt.subplots(r, c, sharex=True)
    for ax, coil in zip(np.atleast_1d(axs).flatten(), pf_active.coil):
        ax.plot(pf_active.time, coil.current.data, color='r')
        tx = ax.twinx()
        tx.plot(pf_active.time, coil.voltage.data, color='b')
        ax.set_title(coil.name)


if __name__ == "__main__":
    import imas
    import sys

    try:
        path = sys.argv[1]
    except IndexError:
        print(f"Usage: {sys.argv[0]} $PATH")
        exit()

    with imas.DBEntry(f"imas:hdf5?path={path}", 'r') as db:
        plot_pf_active(db.get("pf_active"))

    plt.show()
