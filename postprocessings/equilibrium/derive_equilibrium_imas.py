import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.interpolate import RectBivariateSpline as spline2
from scipy.optimize import minimize_scalar
import logging
logger = logging.getLogger(__name__)


def find_rz_contour(
    psi_targets,
    lr,
    lz,
    fpsirz,
    fpsirzr,
    fpsirzz,
    r,
    z,
    maxiter=100,
    under_relax=0.3,
):
    for i in range(maxiter):
        psi = fpsirz(r, z, grid=False)
        if all(np.allclose(psi_, psi_target) for psi_, psi_target in zip(psi, psi_targets)):
            return r, z

        c = (
            under_relax
            * (psi_targets[:, None] - psi)
            / (lr * fpsirzr(r, z, grid=False) + lz * fpsirzz(r, z, grid=False))
        )
        r, z = r + c * lr, z + c * lz

    for psi_, psi_target in zip(psi, psi_targets):
        if not np.allclose(psi_, psi_target):
            logger.warn(f"Contour finding not converged for {psi_target = }")

    raise RuntimeError("Maximum iteration reached")


def find_minmax(r, z):
    fr = make_interp_spline(np.arange(len(r) + 1), np.hstack([r, r[0]]), bc_type="periodic")
    fz = make_interp_spline(np.arange(len(z) + 1), np.hstack([z, z[0]]), bc_type="periodic")
    xrmin = minimize_scalar(fr, bounds=(0, len(r) + 1))
    xrmax = minimize_scalar(lambda x: -fr(x), bounds=(0, len(r) + 1))
    xzmin = minimize_scalar(fz, bounds=(0, len(z) + 1))
    xzmax = minimize_scalar(lambda x: -fz(x), bounds=(0, len(z) + 1))

    return fr(xrmin.x), fr(xrmax.x), fr(xzmin.x), fr(xzmax.x), fz(xrmin.x), fz(xrmax.x), fz(xzmin.x), fz(xzmax.x)


def derive_equilibrium_imas_slice(
    equilibrium,
    it=-1,
    ip2d=0,
    psi_norm_boundary=0.999,
    nth=128,
    maxiter=100,
    under_relax=0.3,
):
    vacuum_toroidal_field = equilibrium.vacuum_toroidal_field
    time_slice = equilibrium.time_slice[it]
    profiles_1d = time_slice.profiles_1d
    profiles_2d = time_slice.profiles_2d[ip2d]
    global_quantities = time_slice.global_quantities
    boundary_separatrix = time_slice.boundary_separatrix

    # Minimal check if psi_axis and psi_boundary is properly defined
    if global_quantities.psi_axis == global_quantities.psi_boundary:
        return

    # Define dpsi and psin, check whether psi contains psi_axis and psi_boundary
    dpsi = global_quantities.psi_boundary - global_quantities.psi_axis
    psin = (profiles_1d.psi - global_quantities.psi_axis)/dpsi
    assert(np.isclose(psin[0], 0.0) and np.isclose(psin[-1], 1.0))

    # Create splines of 2D psi
    rgrid = profiles_2d.grid.dim1
    zgrid = profiles_2d.grid.dim2
    fpsirz = spline2(rgrid, zgrid, (profiles_2d.psi - global_quantities.psi_axis)/dpsi)
    fpsirzr = fpsirz.partial_derivative(dx=1, dy=0)
    fpsirzz = fpsirz.partial_derivative(dx=0, dy=1)

    # Find Br and Bz from partial derivatives of psirz
    profiles_2d.b_field_r = fpsirzz(rgrid, zgrid)*(dpsi/(2.0*np.pi*rgrid[:, None]))
    profiles_2d.b_field_z = -fpsirzr(rgrid, zgrid)*(dpsi/(2.0*np.pi*rgrid[:, None]))

    # Create splines of 1D profile function (only f and q)
    ff = make_interp_spline(psin, profiles_1d.f)
    fq = make_interp_spline(psin, profiles_1d.q)

    # Define raxis, zaxis, baxis, q0, and b0
    raxis = global_quantities.magnetic_axis.r
    zaxis = global_quantities.magnetic_axis.z
    baxis = ff(0.0)/raxis
    q0 = fq(0.0)
    b0 = vacuum_toroidal_field.b0[it]

    # Calculate quantities that does not require flux-averaging (rho_tor_norm is overwritten here)
    profiles_1d.phi = make_interp_spline(psin, profiles_1d.q).antiderivative()(psin)*dpsi
    profiles_1d.rho_tor = (profiles_1d.phi/(np.pi*b0))**0.5
    profiles_1d.rho_tor_norm = profiles_1d.rho_tor/profiles_1d.rho_tor[-1]
    profiles_1d.dpsi_drho_tor = 2.0*np.pi*b0*profiles_1d.rho_tor/profiles_1d.q
    profiles_1d.magnetic_shear = (
        make_interp_spline(psin, np.log(np.abs(profiles_1d.q) + 1e-9)).derivative()(psin)/dpsi*
        profiles_1d.rho_tor*profiles_1d.dpsi_drho_tor
    )

    # Find contour points, except for the magnetic axis
    theta = np.linspace(0.0, 2.0*np.pi, nth, endpoint=False)
    lr = np.cos(theta)
    lz = np.sin(theta)
    try:
        r, z = find_rz_contour(
            np.hstack([psin[1:-1], psi_norm_boundary]),
            lr,
            lz,
            fpsirz,
            fpsirzr,
            fpsirzz,
            raxis + np.outer(profiles_1d.rho_tor[1:], lr),
            zaxis + np.outer(profiles_1d.rho_tor[1:], lz),
            maxiter=maxiter,
            under_relax=under_relax,
        )
    except RuntimeError as e:
        logger.error(f"Failed to find R/Z contour at time = {time_slice.time}")
        return

    # Define 2D quantities
    gradpsi = np.linalg.norm([fpsirzr(r, z, grid=False), fpsirzz(r, z, grid=False)], axis=0)*dpsi
    gradrho = gradpsi/profiles_1d.dpsi_drho_tor[1:, None]
    btor = ff(psin[1:])[:, None]/r
    bpol = gradpsi/(2.0*np.pi*r)
    dl = np.linalg.norm([0.5*(np.roll(r, 1, axis=-1) - np.roll(r, -1, axis=-1)),
                         0.5*(np.roll(z, 1, axis=-1) - np.roll(z, -1, axis=-1))], axis=0)
    jacobian = np.abs(dl/bpol)
    denom = np.sum(jacobian, axis=-1)

    # Calculate min/max values and shaping parameters
    # TODO: Improve b_field calculations
    rrmin, rrmax, rzmin, rzmax, zrmin, zrmax, zzmin, zzmax = np.array([find_minmax(r_, z_) for r_, z_ in zip(r, z)]).T
    bmin = np.linalg.norm(
        [
            ff(psin[1:])/rrmax,
            np.linalg.norm([
                fpsirzr(rrmax, zrmax, grid=False),
                fpsirzz(rrmax, zrmax, grid=False)],
                axis=0
            )*dpsi/(2.0*np.pi*rrmax)
        ],
        axis=0,
    )
    bmax = np.linalg.norm(
        [
            ff(psin[1:])/rrmin,
            np.linalg.norm([
                fpsirzr(rrmin, zrmin, grid=False),
                fpsirzz(rrmin, zrmin, grid=False)],
                axis=0
            )*dpsi/(2.0*np.pi*rrmin)
        ],
        axis=0,
    )
    profiles_1d.r_inboard = np.hstack([raxis, rrmin])
    profiles_1d.r_outboard = np.hstack([raxis, rrmax])
    profiles_1d.b_field_min = np.hstack([np.abs(baxis), bmin])
    profiles_1d.b_field_max = np.hstack([np.abs(baxis), bmax])
    profiles_1d.geometric_axis.r = np.hstack([raxis, (rrmin + rrmax)/2])
    profiles_1d.geometric_axis.z = np.hstack([zaxis, (zzmin + zzmax)/2])
    profiles_1d.elongation = make_interp_spline(psin[1:], (zzmax - zzmin)/(rrmax - rrmin))(psin)
    profiles_1d.triangularity_upper = make_interp_spline(
        psin[1:],
        2.0*((profiles_1d.geometric_axis.r[1:] - rzmax)/(rrmax - rrmin))
    )(psin)
    profiles_1d.triangularity_lower = make_interp_spline(
        psin[1:],
        2.0*((profiles_1d.geometric_axis.r[1:] - rzmin)/(rrmax - rrmin))
    )(psin)

    # Calculate flux-averaged quantities
    profiles_1d.dvolume_dpsi = np.hstack([2.0*np.pi*raxis/baxis*q0, np.sum(dl/bpol, axis=-1)])
    profiles_1d.gm1 = np.hstack([1.0/raxis**2, np.sum(jacobian/r**2, axis=-1)/denom])
    profiles_1d.gm2 = make_interp_spline(psin[1:], np.sum(jacobian*(gradrho/r)**2, axis=-1)/denom)(psin)
    profiles_1d.gm3 = make_interp_spline(psin[1:], np.sum(jacobian*gradrho**2, axis=-1)/denom)(psin)
    profiles_1d.gm4 = np.hstack([1/baxis**2, np.sum(jacobian/(btor**2 + bpol**2), axis=-1)/denom])
    profiles_1d.gm5 = np.hstack([baxis**2, np.sum(jacobian*(btor**2 + bpol**2), axis=-1)/denom])
    profiles_1d.gm6 = make_interp_spline(psin[1:], np.sum(jacobian*gradrho/(btor**2 + bpol**2), axis=-1)/denom)(psin)
    profiles_1d.gm7 = make_interp_spline(psin[1:], np.sum(jacobian*gradrho, axis=-1)/denom)(psin)
    profiles_1d.gm8 = np.hstack([raxis, np.sum(jacobian*r, axis=-1)/denom])
    profiles_1d.gm9 = np.hstack([1.0/raxis, np.sum(jacobian/r, axis=-1)/denom])
    profiles_1d.b_field_average = np.hstack([np.abs(baxis), np.sum(jacobian*(btor**2 + bpol**2)**0.5, axis=-1)/denom])

    # Calculate rest of the part using flux-averaged quantities
    profiles_1d.dvolume_drho_tor = 4.0*np.pi**2*b0*profiles_1d.rho_tor/(profiles_1d.f*profiles_1d.gm1)
    profiles_1d.volume = make_interp_spline(psin, profiles_1d.dvolume_dpsi).antiderivative()(psin)*dpsi
    profiles_1d.rho_volume_norm = profiles_1d.volume/profiles_1d.volume[-1]
    profiles_1d.surface = profiles_1d.dvolume_drho_tor*profiles_1d.gm7
    profiles_1d.darea_dpsi = profiles_1d.gm9*profiles_1d.dvolume_dpsi/(2.0*np.pi)
    profiles_1d.darea_drho_tor = profiles_1d.darea_dpsi*profiles_1d.dpsi_drho_tor
    profiles_1d.area = make_interp_spline(psin, profiles_1d.darea_dpsi).antiderivative()(psin)*dpsi
    profiles_1d.j_tor = -2.0*np.pi/profiles_1d.gm9*(
        profiles_1d.dpressure_dpsi +
        profiles_1d.f_df_dpsi*profiles_1d.gm1/(4.0e-7*np.pi)
    )
    profiles_1d.j_parallel = -2.0*np.pi/b0*(
        profiles_1d.f*profiles_1d.dpressure_dpsi +
        profiles_1d.gm5*profiles_1d.f_df_dpsi/profiles_1d.f/(4.0e-7*np.pi)
    )

    # TODO: squareness, trapped_fraction, beta_pol, mass_density


def derive_equilibrium_imas(equilibrium, *args, **kwargs):
    for it in range(len(equilibrium.time_slice)):
        derive_equilibrium_imas_slice(equilibrium, it, *args, **kwargs)

