import numpy as np
from scipy.interpolate import make_interp_spline
import logging
logger = logging.getLogger(__name__)

try:
    import fit_at_time
    import parameter_settings
except ImportError:
    import os
    import sys
    sys.path.insert(0, os.path.abspath(os.environ["BAYESIAN_FITTING_DIR"]))
    import fit_at_time
    import parameter_settings

BURNIN_NUMBER = 30000
SAMPLING_NUMBER = 50000


def bayesian_fitting_imas(
    core_profiles,
    interferometer,
    thomson_scattering,
    charge_exchange,
    equilibrium,
    time,
    n=129,
    average_time=0.1,
    tci_average_time=0.04,
    **kwargs,
):
    # Resize core_profiles, set time
    core_profiles.time = [time]
    core_profiles.profiles_1d.resize(1)
    core_profiles.profiles_1d[0].time = time

    # Get efit_at_time data from equilibrium IDS
    it_eq = np.argmin(np.abs(equilibrium.time - time))
    efit_at_time = get_efit_at_time_imas(equilibrium, it_eq)

    # Set aliases of equilibrium IDS
    time_slice = equilibrium.time_slice[it_eq]
    profiles_1d = time_slice.profiles_1d
    global_quantities = time_slice.global_quantities
    psi_axis = global_quantities.psi_axis
    psi_boundary = global_quantities.psi_boundary

    # Set normalized psi grid and update core_profiles IDS using equilibrium IDS
    x = np.linspace(0.0, 1.0, n)
    core_profiles.profiles_1d[0].grid.rho_pol_norm = x**0.5
    core_profiles.profiles_1d[0].grid.psi = x*(psi_boundary - psi_axis) + psi_axis
    core_profiles.profiles_1d[0].grid.rho_tor_norm = make_interp_spline(
        (profiles_1d.psi - psi_axis)/(psi_boundary - psi_axis),
        profiles_1d.rho_tor_norm,
    )(x)

    # Get data from interferometer IDS
    mask_tci = (
        (time - 0.5*tci_average_time <= interferometer.time) &
        (interferometer.time < time + 0.5*tci_average_time)
    )

    # Get data from thomson_scattering IDS
    it_ts = np.argmin(np.abs(thomson_scattering.time - time))
    r_ts = np.array([channel.position.r + channel.delta_position.r[it_ts] for channel in thomson_scattering.channel])
    z_ts = np.array([channel.position.z + channel.delta_position.z[it_ts] for channel in thomson_scattering.channel])

    # Get data from charge_exchange IDS
    mask_ces = (
        (time - 0.5*average_time <= charge_exchange.time) &
        (charge_exchange.time < time + 0.5*average_time)
    )
    it_ces = np.argmin(np.abs(charge_exchange.time - time))
    r_ces = np.array([channel.position.r.data[it_ces] for channel in charge_exchange.channel])
    z_ces = np.array([channel.position.z.data[it_ces] for channel in charge_exchange.channel])

    # Fit electron density from interferometer and thomson_scattering IDS and put it into core_profiles IDS
    # Caveat: use only the lower error
    #         no time-averaging for TS
    #         use all TCI channels and they should be stored in a correct order
    n_e_line_average = np.mean(
        [
            channel.n_e_line_average.data[mask_tci] for channel in interferometer.channel
        ], axis=1
    )*1.0e-19
    n_e_line_average_error = np.std(
        [
            channel.n_e_line_average.data[mask_tci] for channel in interferometer.channel
        ], axis=1
    )*1.0e-19
    n_e = np.array([channel.n_e.data[it_ts] for channel in thomson_scattering.channel])*1.0e-19
    n_e_error = np.array([channel.n_e.data_error_lower[it_ts] for channel in thomson_scattering.channel])*1.0e-19
    mean, _ = bayesian_fitting_ne(
        x,
        efit_at_time,
        r_ts,
        z_ts,
        n_e,
        n_e_error,
        n_e_line_average,
        n_e_line_average_error,
        **kwargs,
    )
    core_profiles.profiles_1d[0].electrons.density = mean*1.0e19

    # Fit electron temperature from thomson_scattering IDS and put it into core_profiles IDS
    # Caveat: use only the lower error
    #         no time-averaging for TS
    t_e = np.array([channel.t_e.data[it_ts] for channel in thomson_scattering.channel])*1.0e-3
    t_e_error = np.array([channel.t_e.data_error_lower[it_ts] for channel in thomson_scattering.channel])*1.0e-3
    mean, _ = bayesian_fitting(
        x,
        efit_at_time,
        r_ts,
        z_ts,
        t_e,
        t_e_error,
        kind="te",
        **kwargs,
    )
    core_profiles.profiles_1d[0].electrons.temperature = mean*1.0e3

    # Fit ion temperature from charge_exchange IDS and put it into core_profiles IDS
    # Caveat: use only the lower error
    #         one ion without metadata
    assert(all(len(channel.ion) == 1 for channel in charge_exchange.channel)) # Assume only one ion is measured
    t_i = np.nanmean(
        [
            channel.ion[0].t_i.data[mask_ces] for channel in charge_exchange.channel
        ], axis=1
    )*1.0e-3
    t_i_error = np.nanmean(
        [
            channel.ion[0].t_i.data_error_lower[mask_ces]**2 for channel in charge_exchange.channel
        ], axis=1
    )**0.5*1.0e-3
    mean, _ = bayesian_fitting(
        x,
        efit_at_time,
        r_ces,
        z_ces,
        t_i,
        t_i_error,
        kind="ti",
        **kwargs,
    )
    if not len(core_profiles.profiles_1d[0].ion):
        # TODO: fill element[#].a, z_n, z_ion, and label
        core_profiles.profiles_1d[0].ion.resize(1) 
    core_profiles.profiles_1d[0].ion[0].temperature = mean*1.0e3

    # Fit toroidal velocity from charge_exchange IDS and put it into core_profiles IDS
    # Caveat: use only the lower error
    #         one ion without metadata
    assert(all(len(channel.ion) == 1 for channel in charge_exchange.channel)) # Assume only one ion is measured
    v_t = np.nanmean(
        [
            channel.ion[0].velocity_tor.data[mask_ces] for channel in charge_exchange.channel
        ], axis=1
    )*-1.0e-5
    v_t_error = np.nanmean(
        [
            channel.ion[0].velocity_tor.data_error_lower[mask_ces]**2 for channel in charge_exchange.channel
        ], axis=1
    )**0.5*1.0e-5
    mean, _ = bayesian_fitting(
        x,
        efit_at_time,
        r_ces,
        z_ces,
        v_t,
        v_t_error,
        kind="vt",
        **kwargs,
    )
    if not len(core_profiles.profiles_1d[0].ion):
        # TODO: fill element[#].a, z_n, z_ion, and label
        core_profiles.profiles_1d[0].ion.resize(1) 
    core_profiles.profiles_1d[0].ion[0].velocity.toroidal = mean*-1.0e5


def bayesian_fitting(
    x,
    efit_at_time,
    r,
    z,
    value,
    error,
    kind,
    sol_pos=2.26,
):
    match kind:
        case "te":
            _, result = fit_at_time.fit_te_at_time({
                "current_time": 0.0,
                "efitAtTime": efit_at_time,
                "tsTeAtTime": value,
                "tsTeErrAtTime": error,
                "tsTeUse": (error > 1e-5) | (r > sol_pos),
                "psin_axis": x,
                "params": parameter_settings.te_parameter_settings,
                "burnin_number": BURNIN_NUMBER,
                "sampling_number": SAMPLING_NUMBER,
            })
        case "ti":
            _, result = fit_at_time.fit_ti_at_time({
                "current_time": 0.0,
                "efitAtTime": efit_at_time,
                "cesPos": r,
                "cesTiAtTime": value,
                "cesTiErrAtTime": error,
                "cesTiUse": (error > 1e-5) | (r > sol_pos),
                "psin_axis": x,
                "params": parameter_settings.ti_parameter_settings,
                "burnin_number": BURNIN_NUMBER,
                "sampling_number": SAMPLING_NUMBER,
            })
        case "vt":
            _, result = fit_at_time.fit_vt_at_time({
                "current_time": 0.0,
                "efitAtTime": efit_at_time,
                "cesPos": r,
                "cesVtAtTime": value,
                "cesVtErrAtTime": error,
                "cesVtUse": (error > 1e-5) | (r > sol_pos),
                "psin_axis": x,
                "params": parameter_settings.vt_parameter_settings,
                "burnin_number": BURNIN_NUMBER,
                "sampling_number": SAMPLING_NUMBER,
            })
        case _:
            raise RuntimeError(f"Unidentified bayesian fitting kind: '{kind}'")

    return result[f"{kind}_profile_mean"], result[f"{kind}_profile_std"]


def bayesian_fitting_ne(
    x,
    efit_at_time,
    r,
    z,
    value,
    error,
    tcivalue,
    tcierror,
    sol_pos=2.26,
):
    assert(len(tcivalue) == len(tcierror))

    _, result = fit_at_time.fit_ne_at_time({
        "current_time": 0.0,
        "efitAtTime": efit_at_time,
        "tsNeAtTime": value,
        "tsNeErrAtTime": error,
        "tsNeUse": (error > 1e-5) | (r > sol_pos),
        "tciMean": tcivalue,
        "tciStd": tcierror,
        "tciUse": (tcivalue >= 0.0) & (tcivalue < 12.0),
        "tci_dl": 0.005,
        "psin_axis": x,
        "params": parameter_settings.ne_parameter_settings,
        "burnin_number": BURNIN_NUMBER,
        "sampling_number": SAMPLING_NUMBER,
        "use_fixed_scales": {"core": False, "edge": False}
    })

    return result[f"ne_profile_mean"], result[f"ne_profile_std"]


def get_efit_at_time_imas(equilibrium, it, ip2d=0):
    time_slice = equilibrium.time_slice[it]
    profiles_2d = time_slice.profiles_2d[ip2d]
    psi_axis = time_slice.global_quantities.psi_axis
    psi_boundary = time_slice.global_quantities.psi_boundary

    return {
        "r": profiles_2d.grid.dim1,
        "z": profiles_2d.grid.dim2,
        # "xlim" # doesn't seems to be required
        # "ylim" # doesn't seems to be required
        "rmagx": time_slice.global_quantities.magnetic_axis.r,
        "zmagx": time_slice.global_quantities.magnetic_axis.z,
        "psi": profiles_2d.psi.T,
        "rbdry": time_slice.boundary_separatrix.outline.r,
        "zbdry": time_slice.boundary_separatrix.outline.z,
        "sibdry": psi_boundary,
        "simagx": psi_axis,
        "psin": (profiles_2d.psi.T - psi_axis)/(psi_boundary - psi_axis),
        "raxis": profiles_2d.grid.dim1,
        "zaxis": profiles_2d.grid.dim2,
        **dict(zip(("gridR", "gridZ"), np.meshgrid(profiles_2d.grid.dim1, profiles_2d.grid.dim2))),
    }
