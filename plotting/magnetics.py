import matplotlib.pyplot as plt
import numpy as np


def a2xyz(th, ph, phi):
    return np.dot(
        np.array([[np.cos(phi), -np.sin(phi), 0], [np.sin(phi), np.cos(phi), 0], [0, 0, 1]]),
        np.array([np.cos(th)*np.cos(ph), np.cos(th)*np.sin(ph), -np.sin(th)])
    )


def rphiz2xyz(r, phi, z):
    return np.array([r*np.cos(phi), r*np.sin(phi), z])


def p2xyz(position):
    return rphiz2xyz(position.r, position.phi, position.z)


def plot_2d_flux_loop(flux_loop):
    # Assume flux loop is defined as one toroidal loop
    assert(len(flux_loop.position) == 1)
    plt.plot(flux_loop.position[0].r, flux_loop.position[0].z, 'ys')


def plot_2d_rogowski_coil(rogowski_coil):
    plt.plot(*np.array([[position.r, position.z] for position in rogowski_coil.position]).T)


def plot_2d_b_field_pol_probe(b_field_pol_probe, arrow_length=0.1):
    plt.quiver(
        b_field_pol_probe.position.r,
        b_field_pol_probe.position.z,
        np.cos(b_field_pol_probe.toroidal_angle)*np.cos(b_field_pol_probe.poloidal_angle),
        -np.cos(b_field_pol_probe.toroidal_angle)*np.sin(b_field_pol_probe.poloidal_angle),
        color='k',
    )


def plot_2d(magnetics):
    for flux_loop in magnetics.flux_loop:
        plot_2d_flux_loop(flux_loop)

    for rogowski_coil in magnetics.rogowski_coil:
        plot_2d_rogowski_coil(rogowski_coil)

    for b_field_pol_probe in magnetics.b_field_pol_probe:
        plot_2d_b_field_pol_probe(b_field_pol_probe)


def plot_3d_flux_loop(flux_loop, nphi=65):
    # Assume flux loop is defined as one toroidal loop
    assert(len(flux_loop.position) == 1)
    r, z = flux_loop.position[0].r, flux_loop.position[0].z
    phi = np.linspace(0, 2*np.pi, nphi)
    plt.plot(r*np.cos(phi), r*np.sin(phi), z, 'y')


def plot_3d_rogowski_coil(rogowski_coil):
    r, phi, z = np.array([[position.r, position.phi, position.z] for position in rogowski_coil.position]).T
    plt.plot(*rphiz2xyz(r, phi, z))


def plot_3d_b_field_pol_probe(b_field_pol_probe, arrow_length=0.1):
    xyz = p2xyz(b_field_pol_probe.position)
    uvw = arrow_length*a2xyz(b_field_pol_probe.poloidal_angle, b_field_pol_probe.toroidal_angle, b_field_pol_probe.position.phi)
    plt.quiver(*xyz, *uvw, color='k')


def plot_3d(magnetics):
    for flux_loop in magnetics.flux_loop:
        plot_3d_flux_loop(flux_loop)

    for rogowski_coil in magnetics.rogowski_coil:
        plot_3d_rogowski_coil(rogowski_coil)

    for b_field_pol_probe in magnetics.b_field_pol_probe:
        plot_3d_b_field_pol_probe(b_field_pol_probe)


def plot_magnetics(magnetics):
    n = len(magnetics.flux_loop)
    r = int(np.sqrt(n))
    c = (n + r - 1)//r

    fig, axs = plt.subplots(r, c, sharex=True)
    for ax, flux_loop in zip(np.atleast_1d(axs).flatten(), magnetics.flux_loop):
        ax.plot(magnetics.time, flux_loop.flux.data, color='r')
        ax.set_title(flux_loop.name)

    n = len(magnetics.b_field_pol_probe)
    r = int(np.sqrt(n))
    c = (n + r - 1)//r

    fig, axs = plt.subplots(r, c, sharex=True)
    for ax, b_field_pol_probe in zip(np.atleast_1d(axs).flatten(), magnetics.b_field_pol_probe):
        ax.plot(magnetics.time, b_field_pol_probe.field.data, color='r')
        ax.set_title(b_field_pol_probe.name)

    n = len(magnetics.rogowski_coil)
    r = int(np.sqrt(n))
    c = (n + r - 1)//r

    fig, axs = plt.subplots(r, c, sharex=True)
    for ax, rogowski_coil in zip(np.atleast_1d(axs).flatten(), magnetics.rogowski_coil):
        ax.plot(magnetics.time, rogowski_coil.current.data, color='r')
        ax.set_title(rogowski_coil.name)


if __name__ == "__main__":
    import imas
    import sys

    try:
        path = sys.argv[1]
    except IndexError:
        print(f"Usage: {sys.argv[0]} $PATH")
        exit()

    with imas.DBEntry(f"imas:hdf5?path={path}", 'r') as db:
        plot_magnetics(db.get("magnetics"))

    plt.show()
