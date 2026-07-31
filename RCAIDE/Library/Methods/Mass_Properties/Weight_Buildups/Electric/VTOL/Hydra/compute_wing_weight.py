# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/Hydra/compute_wing_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

# RCAIDE imports
from RCAIDE.Framework.Core import Units
from RCAIDE.Library.Attributes.Materials import AS4_3502_Unidirectional_Carbon_Fiber

# package imports
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  Constants
# ----------------------------------------------------------------------------------------------------------------------
# Direct port of Hydra's `src/Python/Stage_1/wing_class.py::wing_group.weight_estimate` (the
# 'target_freq' branch — the model actually selected by default in the source, `model =
# 'target_freq'`) and its dependencies `motor_mount_mass.py::motor_mount_mass`/`spar_mass`/
# `coefficients_simple`/`bbars`. This is the "class-based" layer `empty_hydra.py` (the SUAVE/
# RCAIDE port audited 2026-07-15) actually descends from — simpler and less tiltrotor-differentiated
# than the `NDARC` method's Chappell-Peyran port (flat 7.5% tilt fraction here vs. `NDARC`'s
# tilt_wing-differentiated 13%). See
# 01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md.
#
# Uses RCAIDE's own `Units` module for SI<->imperial conversion (not Hydra's own rounded
# constants) — validated within tight relative tolerance of the original source, not bit-for-bit
# identical to it, for the same reason documented in `../NDARC/compute_wing_weight.py`.

_LB2KG = Units.lb   # multiply an lb quantity by this to get kg
_KG2LB = 1.0 / Units.lb
_M2F   = 1.0 / Units.ft

_TAU_W  = 0.158   # wing thickness-to-chord ratio
_TAPER  = 0.8     # taper ratio, ctip/croot
_NZ     = 3.8     # max design load factor
_F_MOUNT = 7.0    # motor mount target natural frequency, Hz
_F_LGLOC = 1.0    # landing-gear-on-wing location factor (1.7247 if LG on wing, 1.0 otherwise)

# Wing spar target natural frequency -- MTOW-scaled, not a flat constant. Source Theory Manual
# only gives a qualitative reference point ("at least 4 Hz desirable for a take off mass of
# 2000 kg... larger aircraft feature progressively relaxed natural frequency constraints"), no
# formula. Calibrated 2026-07-28 (00-decisions/2026-07-28-hydra-wing-method-selected-and-
# calibrated.md): fn=2.65 Hz at MTOW=2796.9 kg reproduces CADDEE_alpha's nasa_lpc regression
# (whose own training data -- Ruh et al. 2023, arXiv:2304.14889 -- included boom-mounted rotor
# loads on the wing) on lift_cruise, more trustworthy for this boom-mounted-rotor architecture
# than the source manual's own un-derived 2000 kg reference point. Scaled from there using the
# standard Froude/dynamic-similarity relation, omega ~ mass^(-1/6) (n ~ mass^(1/3) volumetric
# scaling, omega ~ n^(-1/2)).
_FN_REF_HZ   = 2.65
_FN_REF_MTOW = 2796.9  # kg


def _target_frequency(vehicle_mtow):
    """ MTOW-scaled wing spar target natural frequency, Hz. See _FN_REF_HZ/_FN_REF_MTOW above. """
    return _FN_REF_HZ * (_FN_REF_MTOW / vehicle_mtow) ** (1.0 / 6.0)

# Motor-mount / spar structural sizing constants (`motor_mount_mass.py`) -- ORIGINAL Hydra source
# values, now only a last-resort fallback inside `_motor_mount_mass`/`_spar_mass` themselves.
# `compute_wing_weight_group`'s own default (when no `spar_material` is passed) is the real,
# sourced `AS4_3502_Unidirectional_Carbon_Fiber` material below, not these constants -- see
# `_default_spar_material()`, decided 2026-07-30 (same precedent as NDARC's `empirical_2026`
# motor-weight default swap). These constants remain reachable only if a caller manages to get
# `E`/`rho`/`sigma_max`/`tau_max` all the way down to `None` despite `spar_material` never being
# `None` at the `compute_wing_weight_group` level -- effectively dead in normal use.
_E         = 122.0e9   # Young's modulus, Pa
_RHO       = 1650.0    # density, kg/m^3
_SIGMA_MAX = 275.0e6   # unidirectional axial limiting stress, Pa
_TAU_MAX   = 47.0e6    # BID ultimate shear stress, Pa
_G         = 9.81      # m/s^2
_TMIN      = 10e-4     # minimum gauge, m (2 layers of carbon)


def _default_spar_material():
    """ Real, sourced default spar material -- AS4/3502 unidirectional carbon tape, UNNOTCHED
        B-basis (the class's own default fields, 180F/wet). Replaces the original hardcoded
        `_SIGMA_MAX`/`_E`/`_RHO`/`_TAU_MAX` constants as `compute_wing_weight_group`'s actual
        default, 2026-07-30. Deliberately NOT the notched/open-hole value (439 MPa) used
        elsewhere this session -- corrected 2026-07-30, same session: the wing spar (like a rotor
        blade) is modeled here as a continuous laminate beam with no fastener/hole anywhere in
        the sizing formula, unlike the fuselage's landing-bearing bolt-pad term, which is a real
        bolted joint. Open-hole test data is the wrong physical basis for a continuous member's
        allowable. See
        01-mission-profiles/00-decisions/2026-07-30-fuselage-weight-formula-comparison.md. """
    return AS4_3502_Unidirectional_Carbon_Fiber()

# Wing skin / control-surface constants
_TSKIN      = 15e-4   # skin thickness, m (3 layers, 0.5mm each)
_LBYC       = 2.1     # length-to-chord ratio for skin mass estimate
_RHO_SKIN   = 1650.0  # mean skin density, kg/m^3 (folds in foam)
_F_TILT     = 0.075   # tilt-actuator weight fraction of wing weight


def _motor_mount_mass(L, M, r, f, E=None, rho=None):
    """ Mass of a uniform cantilever tube beam with a tip mass, sized to a target natural
        frequency. Direct port of `motor_mount_mass.py::motor_mount_mass`.

        Inputs:
            L   beam length                    [m]
            M   tip mass                       [kg]
            r   outer radius, circular section  [m]
            f   target natural frequency        [Hz]
            E   Young's modulus override, optional -- falls back to module default (`_E`) if None
            rho density override, optional -- falls back to module default (`_RHO`) if None

        Outputs:
            mass  beam mass                    [kg]
    """
    E   = _E   if E   is None else E
    rho = _RHO if rho is None else rho

    wn = f * 2.0 * np.pi

    LHS = 1.0
    fcoef = 11.0 / 420.0
    L4 = L*L*L*L
    A = LHS / (3*L4)

    dRHS = M / (9.0*L)
    coeff = A*E*r*r / (2.0*rho*wn*wn) - fcoef
    mass = dRHS / coeff

    const = 2*np.pi*r*rho

    t = mass / const
    t = max(_TMIN, t)
    mass = t * const
    mass = mass * L
    return mass


def _bbars(taper):
    """ Deflection-term coefficients for a tapered cantilever beam under tip load. Direct port of
        `motor_mount_mass.py::bbars` — a dependency of `_coefficients_simple` only. """
    temp = ((taper+1)/(taper-1))
    temp = temp*temp*temp
    b1bar = temp*(taper-1.0)*(taper-2.0)
    bmbar = taper*temp/(taper+1)
    blbar = temp
    b0bar = -temp*0.5*taper
    brbar = -temp
    return b0bar, brbar, b1bar, blbar, bmbar


def _coefficients_simple(r1, taper, L, masses, xposn, tbyc):
    """ KE/PE expression coefficients for a tapered beam with lumped masses — cubic-polynomial
        curve fits of the exact expressions. Direct port of
        `motor_mount_mass.py::coefficients_simple` — a dependency of `_spar_mass` only. """
    r2 = taper*r1
    rprime = (r2-r1)/L
    ratio = r1/r2  # noqa: F841 (kept for parity with source; not consumed further)

    cbar = r1/tbyc*(1+taper)
    AR = 2*L/cbar
    K = 2.0*tbyc/AR
    K3inv = 1.0/(K*K*K)
    taper2 = taper*taper
    taper3 = taper*taper2

    ABC2 = (1.626090567083019 - 0.482597591190902*taper +
            2.119679981957482*taper2 - 0.595753854491954*taper3)

    LHS = ABC2*K3inv

    logr1 = np.log(r1)

    scale = L*L*K3inv/(K*K)

    RHS_C = (0.369176557170034*taper3 + 0.151779690454288*taper2 +
              0.269000169945435*taper + 0.048225660574691)

    RHS_v3 = RHS_C*scale

    dRHS = 0.0
    b0bar, brbar, b1bar, blbar, bmbar = _bbars(taper)
    b0 = (b0bar+brbar*logr1)*K3inv
    b1 = b1bar*0.5/L*K3inv
    bm = bmbar*L*K*0.5*K3inv
    bl = blbar*K3inv
    for Mk, xk in zip(masses, xposn):
        r = r1 + rprime*xk
        w = b0 + b1*xk + bm/r + bl*np.log(r)
        dRHS = dRHS + Mk*0.5*w*w

    return LHS, RHS_v3, dRHS


def _spar_mass(r1, taper, L, m_by_a, masses, xposn, thrusts, fn=6.0, tbyc=0.168,
               E=None, rho=None, sigma_max=None, tau_max=None):
    """ Mass of a tapered cantilever spar sized to a target natural frequency, with distributed
        (skin) and lumped (motor mount/motor/rotor) non-structural mass, checked against static
        stress limits in hover. Direct port of `motor_mount_mass.py::spar_mass`.

        Inputs:
            r1        root cross-section tube radius            [m]
            taper     taper ratio, tip chord/root chord         [-]
            L         beam length (half-span)                   [m]
            m_by_a    non-structural mass per unit area          [kg/m^2]
            masses    lumped masses (mounts, motors, rotors)     [kg]
            xposn     spanwise distance of each lumped mass      [m]
            thrusts   rotor thrust at each lumped mass position  [N]
            fn        target natural frequency, half-wing        [Hz]
            tbyc      wing thickness-to-chord ratio              [-]
            E         Young's modulus override, optional -- falls back to `_E` if None
            rho       density override, optional -- falls back to `_RHO` if None
            sigma_max axial stress allowable override, optional -- falls back to `_SIGMA_MAX`
            tau_max   shear stress allowable override, optional -- falls back to `_TAU_MAX`

        Outputs:
            Mspar     spar mass                                 [kg]
    """
    E         = _E         if E         is None else E
    rho       = _RHO       if rho       is None else rho
    sigma_max = _SIGMA_MAX if sigma_max is None else sigma_max
    tau_max   = _TAU_MAX   if tau_max   is None else tau_max

    LHS, RHS, dRHS = _coefficients_simple(r1, taper, L, masses, xposn, tbyc)
    omegan = fn*2*np.pi
    beta = m_by_a/tbyc
    coef = E*np.pi*0.5/(omegan*omegan)*LHS - rho*np.pi*RHS
    b = beta*RHS + dRHS
    t = b/coef
    rbar = (1+taper)*0.5*r1

    t = max(t, _TMIN)

    xnodes = list(xposn)
    lumped_masses = list(masses)
    lumped_thrusts = list(thrusts)
    xnodes.insert(0, 0.0)
    lumped_masses.insert(0, 0.0)
    lumped_thrusts.insert(0, 0.0)
    nelem = len(xnodes)-1
    nsample = 5
    xpts = np.zeros((nsample*nelem))
    ipt = 0
    BM = np.zeros_like(xpts)
    Shear = np.zeros_like(xpts)
    for i in range(nelem):
        dx = (xnodes[i+1] - xnodes[i])/float(nsample)
        for j in range(nsample):
            xpts[ipt] = xnodes[i] + dx*float(j)

            for x, M, T in zip(xnodes, lumped_masses, lumped_thrusts):
                if x > xpts[ipt]:
                    Fz = T - M*_G
                    BM[ipt] = BM[ipt] + Fz*(x - xpts[ipt])
                    Shear[ipt] = Shear[ipt] + Fz
            ipt = ipt + 1

    rprime = r1*(taper-1.0)/L
    radius = r1 + rprime*xpts

    sigma = np.divide(BM*_NZ, 8*np.pi*t*np.square(radius))
    SF = np.amin(sigma_max/sigma)

    tshear = t
    CSArea = 2*np.pi*radius*tshear
    tau = 2.0*Shear/CSArea*_NZ
    SF2 = np.amin(tau_max/tau)

    if SF < 1.5:
        t = t*1.5/SF
    if SF2 < 1.5:
        tshear = tshear*1.5/SF2

    mspar = 2*np.pi*rbar*(t+tshear)*rho
    Mspar = mspar*L

    return Mspar


def compute_wing_weight_group(vehicle_mtow, n_wings, aspect_ratio, area, lift_fraction,
                               rotor_radius, rotor_mass_assembly, rotor_y_positions,
                               rotor_masses, rotor_thrusts, wing_flap_redundancy=1.0,
                               tilt_actuator_redundancy=1.0, tech_factor_wing=1.0,
                               tech_factor_flight_control=1.0, spar_material=None):
    """ Calculates the structural mass of one fixed-wing group using Hydra's target-frequency
        spar-sizing method (the class-based layer `empty_hydra.py` descends from).

        Sizes the motor mounts (cantilever tube, tip mass, target frequency) and the wing spar
        itself (tapered cantilever beam, target frequency, distributed skin mass + lumped
        mount/motor/rotor masses, checked against static stress in hover) independently, then adds
        control-surface actuator weight (scaling law) and a flat tilt-actuator fraction of wing
        weight — no `tilt_wing`/`tilt_rotor` distinction (unlike the `NDARC` method's
        Chappell-Peyran port, this layer's tilt term is a single flat fraction regardless of
        mechanism type; confirmed against source 2026-07-15 audit).

        Inputs:
            vehicle_mtow                  vehicle max takeoff weight                    [kg]
            n_wings                       number of wings in this group (e.g. 2 biplane) [Unitless]
            aspect_ratio                  wing group aspect ratio                       [Unitless]
            area                          wing group reference area, one wing           [m^2]
            lift_fraction                 fraction of MTOW-derived lift this group carries, 0-1
            rotor_radius                  radius of the rotor type on this wing         [m]
            rotor_mass_assembly           mass of one rotor+hub+actuator assembly       [kg]
            rotor_y_positions             spanwise location of each rotor, per half-wing [m]
            rotor_masses                  motor+rotor mass at each location (no mount)   [kg]
            rotor_thrusts                 thrust at each location                        [N]
            wing_flap_redundancy          actuator weight redundancy multiplier          [Unitless]
            tilt_actuator_redundancy      tilt actuator weight redundancy multiplier     [Unitless]
            tech_factor_wing              technology factor, wing structure              [Unitless]
            tech_factor_flight_control    technology factor, actuators/tilters           [Unitless]
            spar_material                 optional `Solid` material instance (e.g.
                                            `AS4_3502_Unidirectional_Carbon_Fiber`). Defaults
                                            (2026-07-30) to `_default_spar_material()` -- a real,
                                            sourced AS4/3502 instance, notched basis -- NOT the
                                            original hardcoded `_E`/`_RHO`/`_SIGMA_MAX`/`_TAU_MAX`
                                            constants (this is a real default-behavior change, not
                                            backward-compatible; those constants are now only a
                                            last-resort fallback deep inside `_motor_mount_mass`/
                                            `_spar_mass`, effectively unreachable in normal use).
                                            Reads `.ultimate_tensile_strength` (sigma_max),
                                            `.ultimate_shear_strength` (tau_max), `.density` (rho),
                                            and the non-standard `.youngs_modulus` attribute if
                                            present (`Solid`'s base interface has no modulus
                                            field).

        Outputs:
            weight:   dict with 'structure', 'actuators', 'tilters', 'mounts', all         [kg]
    """
    if spar_material is None:
        spar_material = _default_spar_material()
    mat_E   = getattr(spar_material, 'youngs_modulus',        None)
    mat_rho = getattr(spar_material, 'density',                None)
    mat_sig = getattr(spar_material, 'ultimate_tensile_strength', None)
    mat_tau = getattr(spar_material, 'ultimate_shear_strength',   None)

    span = np.sqrt(aspect_ratio * area)
    chord = area / span

    W = vehicle_mtow * _KG2LB
    Sw = area * _M2F * _M2F
    fL = lift_fraction
    Wt_wing = fL * W / n_wings  # noqa: F841 (kept for parity with source; not consumed further)

    nrotors_half_wing = len(rotor_y_positions)

    L_mount = rotor_radius + chord*0.3
    r_tube = 0.125*rotor_radius*0.5
    mount_mass = _motor_mount_mass(L_mount, rotor_mass_assembly, r_tube, _F_MOUNT,
                                    E=mat_E, rho=mat_rho)

    mounts_wt = _KG2LB * mount_mass * nrotors_half_wing * (2*n_wings)

    m_by_a = _TSKIN * _LBYC * _RHO_SKIN

    rbar = _TAU_W*chord*0.5
    rroot = 2*rbar/(1+_TAPER)
    L = span*0.5

    y = np.asarray(rotor_y_positions)
    Mk = np.asarray(rotor_masses) + mount_mass
    T = np.asarray(rotor_thrusts)
    M_spar = _spar_mass(rroot, _TAPER, L, m_by_a, Mk, y, T, _target_frequency(vehicle_mtow), _TAU_W,
                         E=mat_E, rho=mat_rho, sigma_max=mat_sig, tau_max=mat_tau)
    # NOTE: `area` (m^2) is used directly here, not converted to ft^2 — matches the original
    # source exactly (`MbyA*self.area`, not `MbyA*Sw`); the trailing *2.2 converts the resulting
    # kg subtotal to lb (a hardcoded approximate kg->lb factor in the source, not the module's own
    # exact constant — kept literal for fidelity).
    wt = (m_by_a*area + M_spar*2)*2.2

    wt = wt * n_wings
    wing_wt = wt
    area_total = Sw*n_wings

    actuator_wt = 0.01735*(W**0.6435)*(area_total**0.40952)
    actuator_wt = actuator_wt * wing_flap_redundancy

    tilt_wt = _F_TILT*wing_wt*tilt_actuator_redundancy

    wing_wt = wing_wt * tech_factor_wing
    actuator_wt = actuator_wt * tech_factor_flight_control
    tilt_wt = tilt_wt * tech_factor_flight_control

    return {
        'structure': wing_wt * _LB2KG,
        'actuators': actuator_wt * _LB2KG,
        'tilters': tilt_wt * _LB2KG,
        'mounts': mounts_wt * _LB2KG,
    }
