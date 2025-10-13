import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np


def plot_ece(ece, time):
    idx = np.argmin(np.abs(ece.time - time))

    fig = plt.figure()
    gs = gridspec.GridSpec(1, 3)
    ax1 = fig.add_subplot(gs[0, :2])
    ax2 = fig.add_subplot(gs[0, 2])

    for channel in ece.channel:
        ax1.plot(ece.time, channel.t_e.data)
    ax1.axvline(ece.time[idx], color='k', linestyle='--')
    ax1.set_xlabel(r"Time [$s$]")
    ax1.set_xlabel(r"Electron temperature [$eV$]")

    r = [channel.position.r[idx] for channel in ece.channel]
    t_e = [channel.t_e.data[idx] for channel in ece.channel]
    ax2.scatter(r, t_e, marker='x')
    ax2.set_xlabel(r"R [$m$]")
    ax2.set_ylabel(r"Electron temperature [$eV$]")


if __name__ == "__main__":
    import imas
    import sys

    try:
        path = sys.argv[1]
        time = float(sys.argv[2])
    except IndexError:
        print(f"Usage: {sys.argv[0]} $PATH $TIME")
        exit()

    with imas.DBEntry(f"imas:hdf5?path={path}", 'r') as db:
        plot_ece(db.get("ece"), time)

    plt.show()
