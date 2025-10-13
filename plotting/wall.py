import matplotlib.pyplot as plt
import numpy as np


def plot_2d_description_2d(description_2d):
    for unit in description_2d.limiter.unit:
        r, z = unit.outline.r, unit.outline.z
        plt.plot(r, z, 'kx-')

    if description_2d.vessel.type.index == 0:
        for unit in description_2d.vessel.unit:
            r, z = unit.annular.centreline.r, unit.annular.centreline.z
            plt.plot(r, z, 'gray')


def plot_2d(wall):
    for description_2d in wall.description_2d:
        plot_2d_description_2d(description_2d)


def plot_3d_description_2d(description_2d, nphi=65):
    phi = np.linspace(0, 2*np.pi, nphi)

    for unit in description_2d.limiter.unit:
        r, z = unit.outline.r, unit.outline.z
        plt.gca().plot_surface(
            np.outer(r, np.cos(phi)),
            np.outer(r, np.sin(phi)),
            np.array([z]*nphi).T,
            color="k",
            alpha=0.2,
        )

    if description_2d.vessel.type.index == 0:
        for unit in description_2d.vessel.unit:
            r, z = unit.annular.centreline.r, unit.annular.centreline.z
            plt.gca().plot_surface(
                np.outer(r, np.cos(phi)),
                np.outer(r, np.sin(phi)),
                np.array([z]*nphi).T,
                color="gray",
                alpha=0.1,
            )


def plot_3d(wall):
    for description_2d in wall.description_2d:
        plot_3d_description_2d(description_2d)
