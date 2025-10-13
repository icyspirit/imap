import matplotlib.pyplot as plt
import numpy as np


def rzphi2xyz(rzphi):
    return np.array([rzphi.r[0]*np.cos(rzphi.phi[0]), rzphi.r[0]*np.sin(rzphi.phi[0]), rzphi.z[0]])


def plot_2d_beam(beam, color=None, l=1.5, n=65):
    th, phi = beam.steering_angle_pol[0], beam.steering_angle_tor[0]
    k = np.array([
        -np.cos(th)*np.cos(phi),
        np.sin(phi),
        -np.sin(th)*np.cos(phi)
    ])
    c = np.array([
        k[0]*np.cos(beam.launching_position.phi[0]) - k[1]*np.sin(beam.launching_position.phi[0]),
        k[0]*np.sin(beam.launching_position.phi[0]) + k[1]*np.cos(beam.launching_position.phi[0]),
        k[2]
    ])
    launching_position = rzphi2xyz(beam.launching_position)
    trajectory = launching_position[:, None] + np.outer(c, np.linspace(0, l, n))
    plt.plot(np.linalg.norm(trajectory[:2], axis=0), trajectory[2], color=color)


def plot_2d(ec_launchers):
    for beam, color in zip(ec_launchers.beam, ['r', 'g', 'b', 'm', 'k', 'c']):
        plot_2d_beam(beam, color)


def plot_3d_beam(beam, color=None, l=1.5):
    th, phi = beam.steering_angle_pol[0], beam.steering_angle_tor[0]
    k = np.array([
        -np.cos(th)*np.cos(phi),
        np.sin(phi),
        -np.sin(th)*np.cos(phi)
    ])
    c = np.array([
        k[0]*np.cos(beam.launching_position.phi[0]) - k[1]*np.sin(beam.launching_position.phi[0]),
        k[0]*np.sin(beam.launching_position.phi[0]) + k[1]*np.cos(beam.launching_position.phi[0]),
        k[2]
    ])
    launching_position = rzphi2xyz(beam.launching_position)
    plt.plot(*np.vstack([launching_position, launching_position + c*l]).T, color=color)


def plot_3d(ec_launchers):
    for beam, color in zip(ec_launchers.beam, ['r', 'g', 'b', 'm', 'k', 'c']):
        plot_3d_beam(beam, color)


def plot_ec_launchers(ec_launchers):
    n = len(ec_launchers.beam)
    r = int(np.sqrt(n))
    c = (n + r - 1)//r

    fig, axs = plt.subplots(r, c, sharex=True)
    for ax, beam in zip(np.atleast_1d(axs).flatten(), ec_launchers.beam):
        ax.set_ylabel(r"Power [$W$]")
        ax.plot(ec_launchers.time, beam.power_launched.data, color='r')

        tx = ax.twinx()
        tx.set_ylabel(r"Frequency [$Hz$]")
        if len(beam.frequency.data) > 1:
            tx.plot(ec_launchers.time, beam.frequency.data, color='b')
        else:
            tx.axhline(beam.frequency.data[0], color='b')
        ax.set_title(beam.identifier)


if __name__ == "__main__":
    import imas
    import sys

    try:
        path = sys.argv[1]
    except IndexError:
        print(f"Usage: {sys.argv[0]} $PATH")
        exit()

    with imas.DBEntry(f"imas:hdf5?path={path}", 'r') as db:
        plot_ec_launchers(db.get("ec_launchers"))

    plt.show()
