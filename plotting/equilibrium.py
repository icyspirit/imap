import matplotlib.pyplot as plt
import numpy as np


markevery = 11
axs11 = None
axs12 = None
axs13 = None
axs14 = None
axs15 = None
axs16 = None
axs17 = None
fig2 = None
axs2 = None

def get(obj, path):
    s = path.split('.')
    if len(s) == 1:
        return getattr(obj, path)
    return get(getattr(obj, s[0]), '.'.join(s[1:]))


def print_equilibrium(equilibrium, time):
    tidx = np.argmin(np.abs(equilibrium.time - time))

    for item in [
        "beta_pol",
        "beta_tor",
        "beta_normal",
        "ip",
        "li_3",
        "volume",
        "area",
        "surface",
        "length_pol",
        "psi_axis",
        "psi_boundary",
        "magnetic_axis.r",
        "magnetic_axis.z",
        "magnetic_axis.b_field_tor",
        "current_centre.r",
        "current_centre.z",
        "current_centre.velocity_z",
        "q_axis",
        "q_95",
        "q_min.value",
        "q_min.rho_tor_norm",
        "energy_mhd",
        "psi_external_average",
        "v_external",
        "plasma_inductance",
        "plasma_resistance",
    ]:
        print(item, get(equilibrium.time_slice[tidx].global_quantities, item))


def plot_profiles_1d(profiles_1d):
    global markevery
    markevery = markevery + 2

    x = (profiles_1d.psi - profiles_1d.psi[0])/(profiles_1d.psi[-1] - profiles_1d.psi[0])

    def plot_(axs, i, item):
        axs[i].set_title(item)
        if len(get(profiles_1d, item)):
            axs[i].plot(x, get(profiles_1d, item), 'o-', markevery=markevery)

    def plot(axs, items):
        for i, item in enumerate(items):
            plot_(axs, i, item)

    plot(axs11.flatten(), ["psi", "phi", "pressure", "f", "dpressure_dpsi", "f_df_dpsi"])
    plot(axs12.flatten(), ["j_tor", "j_parallel", "q", "magnetic_shear", "r_inboard", "r_outboard"])
    plot(axs13.flatten(), ["rho_tor", "rho_tor_norm", "dpsi_drho_tor", "geometric_axis.r", "geometric_axis.z"])
    plot(axs14.flatten(), ["elongation", "triangularity_upper", "triangularity_lower", "squareness_upper_inner", "squareness_upper_outer", "squareness_lower_inner", "squareness_lower_outer"])
    plot(axs15.flatten(), ["volume", "rho_volume_norm", "dvolume_dpsi", "area", "darea_dpsi", "darea_drho_tor", "surface", "trapped_fraction"])
    plot(axs16.flatten(), ["gm1", "gm2", "gm3", "gm4", "gm5", "gm6", "gm7", "gm8", "gm9"])
    plot(axs17.flatten(), ["b_field_average", "b_field_min", "b_field_max", "beta_pol", "mass_density"])


def plot_profiles_2d(profiles_2d):
    r, z = np.meshgrid(profiles_2d.grid.dim1, profiles_2d.grid.dim2)

    if profiles_2d.psi.size > 0:
        c = axs2[0].contour(r, z, profiles_2d.psi.T)
        axs2[0].axis("scaled")
        fig2.colorbar(c, ax=axs2[0])

    if profiles_2d.b_field_r.size > 0 and profiles_2d.b_field_z.size > 0:
        axs2[1].quiver(r, z, profiles_2d.b_field_r.T, profiles_2d.b_field_z.T, scale=5)
        axs2[1].axis("scaled")


def plot_equilibrium(equilibrium, time):
    tidx = np.argmin(np.abs(equilibrium.time - time))

    plot_profiles_1d(equilibrium.time_slice[tidx].profiles_1d)
    # Plot only the first index of the profiles_2d (profiles_2d[0])
    if equilibrium.time_slice[tidx].profiles_2d:
        plot_profiles_2d(equilibrium.time_slice[tidx].profiles_2d[0])


if __name__ == "__main__":
    import imas
    import sys

    figsize = (12, 6)
    _, axs11 = plt.subplots(2, 3, sharex=True, figsize=figsize)
    _, axs12 = plt.subplots(2, 3, sharex=True, figsize=figsize)
    _, axs13 = plt.subplots(2, 3, sharex=True, figsize=figsize)
    _, axs14 = plt.subplots(2, 4, sharex=True, figsize=figsize)
    _, axs15 = plt.subplots(2, 4, sharex=True, figsize=figsize)
    _, axs16 = plt.subplots(3, 3, sharex=True, figsize=figsize)
    _, axs17 = plt.subplots(2, 3, sharex=True, figsize=figsize)
    fig2, axs2 = plt.subplots(1, 2, figsize=figsize)

    try:
        paths = sys.argv[1:-1]
        time = float(sys.argv[-1])
    except ValueError:
        print(f"Usage: {sys.argv[0]} $PATHS... $TIME")
        exit()

    for path in paths:
        with imas.DBEntry(f"imas:hdf5?path={path}", 'r') as db:
            equilibrium = db.get("equilibrium")

        print_equilibrium(equilibrium, time)
        plot_equilibrium(equilibrium, time)

    plt.show()

