import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np


def plot_2d_channel(channel):
    # Plotting only the first r, z position data
    if len(channel.position.r.data) and len(channel.position.z.data):
        plt.scatter(channel.position.r.data[0], channel.position.z.data[0], color='r', marker='.')


def plot_2d(charge_exchange):
    for channel in charge_exchange.channel:
        plot_2d_channel(channel)


def plot_charge_exchange(charge_exchange, time):
    if charge_exchange.time.size == 0:
        print("No data")
        return

    idx = np.argmin(np.abs(charge_exchange.time - time))

    fig = plt.figure()
    gs = gridspec.GridSpec(3, 2)
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, :])
    ax3 = fig.add_subplot(gs[2, 0])
    ax4 = fig.add_subplot(gs[2, 1])

    for channel in charge_exchange.channel:
        for ion in channel.ion:
            ax1.plot(charge_exchange.time, ion.t_i.data)
            ax2.plot(charge_exchange.time, ion.velocity_tor.data)
    ax1.axvline(charge_exchange.time[idx], color='k', linestyle='--')
    ax2.axvline(charge_exchange.time[idx], color='k', linestyle='--')
    ax1.set_xlabel(r"Time [$s$]")
    ax2.set_xlabel(r"Time [$s$]")
    ax1.set_ylabel(r"Ion temperature [$eV$]")
    ax2.set_ylabel(r"Toroidal velocity [$m/s$]")

    r = [channel.position.r.data[idx] for channel in charge_exchange.channel]
    t_i = np.array([[ion.t_i.data[idx] for ion in channel.ion] for channel in charge_exchange.channel])
    t_i_upper = np.array([[ion.t_i.data_error_upper[idx] for ion in channel.ion] for channel in charge_exchange.channel])
    t_i_lower = np.array([[ion.t_i.data_error_lower[idx] for ion in channel.ion] for channel in charge_exchange.channel])
    velocity_tor = np.array([[ion.velocity_tor.data[idx] for ion in channel.ion] for channel in charge_exchange.channel])
    velocity_tor_upper = np.array([[ion.velocity_tor.data_error_upper[idx] for ion in channel.ion] for channel in charge_exchange.channel])
    velocity_tor_lower = np.array([[ion.velocity_tor.data_error_lower[idx] for ion in channel.ion] for channel in charge_exchange.channel])

    for v, l, u in zip(t_i.T, t_i_lower.T, t_i_upper.T):
        ax3.errorbar(r, v, yerr=[l, u])
    for v, l, u in zip(velocity_tor.T, velocity_tor_lower.T, velocity_tor_upper.T):
        ax4.errorbar(r, v, yerr=[l, u])
    ax3.set_xlabel(r"R [$m$]")
    ax4.set_xlabel(r"R [$m$]")
    ax3.set_ylabel(r"Ion temperature [$eV$]")
    ax4.set_ylabel(r"Toroidal velocity [$m/s$]")


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
        plot_charge_exchange(db.get("charge_exchange"), time)

    plt.show()
