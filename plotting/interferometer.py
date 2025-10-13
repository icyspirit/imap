import matplotlib.pyplot as plt
import numpy as np


def plot_interferometer(interferometer):
    n = len(interferometer.channel)
    r = int(np.sqrt(n))
    c = (n + r - 1)//r

    fig, axs = plt.subplots(r, c, sharex=True)
    for ax, channel in zip(np.atleast_1d(axs).flatten(), interferometer.channel):
        ax.set_title(channel.name)
        ax.plot(interferometer.time, channel.n_e_line_average.data)


if __name__ == "__main__":
    import imas
    import sys

    try:
        path = sys.argv[1]
    except IndexError:
        print(f"Usage: {sys.argv[0]} $PATH")
        exit()

    with imas.DBEntry(f"imas:hdf5?path={path}", 'r') as db:
        plot_interferometer(db.get("interferometer"))

    plt.show()
