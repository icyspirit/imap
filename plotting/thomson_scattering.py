import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np


def plot_2d_channel(channel):
    plt.scatter(channel.position.r, channel.position.z, color='b', marker='.')


def plot_2d(thomson_scattering):
    for channel in thomson_scattering.channel:
        plot_2d_channel(channel)


def plot_thomson_scattering(thomson_scattering, time):
    idx = np.argmin(np.abs(thomson_scattering.time - time))

    fig = plt.figure()
    gs = gridspec.GridSpec(3, 2)
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, :])
    ax3 = fig.add_subplot(gs[2, 0])
    ax4 = fig.add_subplot(gs[2, 1])

    for channel in thomson_scattering.channel:
        ax1.plot(thomson_scattering.time, channel.n_e.data)
        ax2.plot(thomson_scattering.time, channel.t_e.data)
    ax1.axvline(thomson_scattering.time[idx], color='k', linestyle='--')
    ax2.axvline(thomson_scattering.time[idx], color='k', linestyle='--')
    ax1.set_xlabel(r"Time [$s$]")
    ax2.set_xlabel(r"Time [$s$]")
    ax1.set_ylabel(r"Electron density [$m^{-3}$]")
    ax2.set_ylabel(r"Electron temperature [$eV$]")

    r = [channel.position.r + channel.delta_position.r[idx] for channel in thomson_scattering.channel]
    n_e = [channel.n_e.data[idx] for channel in thomson_scattering.channel]
    n_e_upper = [channel.n_e.data_error_upper[idx] for channel in thomson_scattering.channel]
    n_e_lower = [channel.n_e.data_error_lower[idx] for channel in thomson_scattering.channel]
    t_e = [channel.t_e.data[idx] for channel in thomson_scattering.channel]
    t_e_upper = [channel.t_e.data_error_upper[idx] for channel in thomson_scattering.channel]
    t_e_lower = [channel.t_e.data_error_lower[idx] for channel in thomson_scattering.channel]

    ax3.errorbar(r, n_e, yerr=[n_e_lower, n_e_upper])
    ax4.errorbar(r, t_e, yerr=[t_e_lower, t_e_upper])
    ax3.set_xlabel(r"R [$m$]")
    ax4.set_xlabel(r"R [$m$]")
    ax3.set_ylabel(r"Electron density [$m^{-3}$]")
    ax4.set_ylabel(r"Electron temperature [$eV$]")


if __name__ == "__main__":
    import imas
    import sys

    try:
        path = sys.argv[1]
        time = float(sys.argv[2])
    except (IndexError, ValueError):
        print(f"Usage: {sys.argv[0]} $PATH $TIME")
        exit()

    with imas.DBEntry(f"imas:hdf5?path={path}", 'r') as db:
        plot_thomson_scattering(db.get("thomson_scattering"), time)

    plt.show()
