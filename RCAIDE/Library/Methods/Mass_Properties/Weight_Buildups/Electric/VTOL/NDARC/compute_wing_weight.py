# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/NDARC/compute_wing_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

# RCAIDE imports
from RCAIDE.Framework.Core import Units

# package imports
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  Constants — Section 29-1.1, NDARC Theory Manual v1.11 pg. 245, Chappell & Peyran
# ----------------------------------------------------------------------------------------------------------------------
# Kept in the original imperial units the regression was calibrated in (ft, lb, slug); SI inputs
# are converted internally, via RCAIDE's own `Units` module, so the RCAIDE-facing signature stays
# SI. Direct port of Hydra's `src/Python/Stage_1/afdd/wing.py::tiltrotor_wing_wt` — validated
# against the original source (999-software/hydra/) across tiltrotor/tiltwing branches, multi-wing
# groups, and non-unity tech factors. Note: Hydra's own unit conversion constants
# (`afdd/conversions.py`) are rounded (lb2kg=0.4536) rather than exact (0.45359237) — using
# RCAIDE's exact `Units.lb` here means output is within ~2e-5 relative tolerance of Hydra's own
# rounded-constant result, not bit-for-bit identical to it. See
# 01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md.

_M2F = 1.0 / Units.ft   # meters -> feet
_KG2LB = 1.0 / Units.lb # kg -> lb
_GRAV_FTS2 = 32.2

_G_TB          = 4e6 * 144.0    # torque box shear modulus, lb/sq.ft
_E_TB          = 10e6 * 144.0   # torque box Young's modulus, lb/sq.ft
_E_SP          = _E_TB          # wing spar modulus, lb/sq.ft
_DENS_TB       = 3.1045         # torsion box density, slug/ft3 (carbon fiber)
_DENS_SPAR     = 3.1045         # spar density, slug/ft3 (carbon fiber)
_E_TB_EFF      = 0.8            # structural efficiency factor, torsion box
_E_SP_EFF      = 0.8            # structural efficiency factor, spar
_CT            = 0.75           # weight correction for spar taper (equivalent stiffness)
_CHORD_TB_WING = 0.6            # ratio of torque box chord to wing chord
_OMEGA_B       = 0.5            # wing flap bending freq, normalized by rotor speed, /rev
_OMEGA_C       = 0.8            # wing chordwise bending freq, normalized by rotor speed, /rev
_OMEGA_T       = 0.9            # wing torsion freq, normalized by rotor speed, /rev
_WING_TC       = 0.23           # wing thickness-to-chord ratio

_LOADING_FAIR = 2.0   # unit weight of LE/TE fairings, lb/sq.ft
_LOADING_FLAP = 3.0   # unit weight of control surface, lb/sq.ft
_F_FIT        = 0.12  # weight fraction of the wing (fittings)
_F_FOLD       = 0.0   # wing fold fraction — no fold mechanism modelled (fixed at 0 in source too)
_F_TILT       = 0.13  # tilt-mechanism weight fraction (same coefficient for tiltwing and tiltrotor)


def compute_wing_weight_group(span, area, chord, rotor_angular_velocity, rotor_radius,
                               wing_tip_mass, vehicle_mtow, tilt_wing, tech_factor=1.0,
                               n_wings=1):
    """ Calculates the structural mass of one wing group of a tiltrotor/tiltwing eVTOL using the
        Chappell-Peyran method (NDARC Theory Manual v1.11, Section 29-1.1).

        Unlike RCAIDE's native `Physics_Based` beam-bending method, this sizes the wing from an
        assumed torsion/bending natural-frequency target driven by the tip mass (rotor + nacelle +
        motor group), not from a lift-distribution load case — the appropriate method for a wing
        whose dominant structural driver is a heavy rotating mass at (or near) the tip, and the
        only one of RCAIDE's weight methods that differentiates tiltwing from tiltrotor tilt
        mechanisms.

        Assumptions:
            Torque box and spar both carbon fiber (0/90 layup assumed for torque box modulus).
            No wing fold mechanism (f_fold fixed at 0, matching the Hydra source).

        Source:
            NDARC Theory Manual v1.11, Section 29-1.1 (Chappell & Peyran).
            Ported from `afdd/wing.py::tiltrotor_wing_wt`,
            `github.com/VahanaOpenSource/vtol_sizing`.

        Inputs:
            span                        wing group span                                    [m]
            area                        wing group reference area                          [m^2]
            chord                       wing group mean chord                              [m]
            rotor_angular_velocity      angular velocity of the rotor driving this wing's
                                         tip mass                                           [rad/s]
            rotor_radius                radius of the rotor driving this wing's tip mass    [m]
            wing_tip_mass               mass of the tilting nacelle/rotor/motor group at
                                         this wing's tip                                    [kg]
            vehicle_mtow                vehicle max takeoff weight                          [kg]
            tilt_wing                   True: tiltwing (whole wing tilts).
                                         False: tiltrotor (only the wingtip nacelle tilts)   [bool]
            tech_factor                 technology weight-scaling factor                    [Unitless]
            n_wings                     number of wings in this group (e.g. 2 for biplane)   [Unitless]

        Outputs:
            weight:                     Data dictionary with 'primary', 'fairing', 'flaps',
                                         'fitting', 'tilt', 'folding', 'total'                [kg]
    """
    b_w = span * _M2F
    S_w = area * _M2F * _M2F
    c   = chord * _M2F
    R   = rotor_radius * _M2F
    Omega = rotor_angular_velocity

    W    = vehicle_mtow * _KG2LB
    wtip = wing_tip_mass * _KG2LB

    S_fair  = 0.10 * S_w
    S_flap  = 0.25 * S_w
    r_pylon = 0.3 * R
    mtip    = wtip / _GRAV_FTS2

    f_tip  = wtip / W
    f_mode = 1.0 - f_tip

    t_w  = _WING_TC * c
    c_tb = _CHORD_TB_WING * c

    wtb  = _CHORD_TB_WING
    tauw = _WING_TC

    F_B = (0.073 * np.sin(2*np.pi*(tauw-0.151)/0.1365) + 0.14598*tauw +
           0.610 * np.sin(2*np.pi*(wtb+0.08)/2.1560) - (0.4126-1.6309*tauw) *
           (wtb-0.131) + 0.0081)

    F_C = 0.640424*wtb*wtb - 0.89717*wtb + 0.4615*tauw + 0.655317

    F_T = (((0.27-tauw)/0.12)*0.12739*(-0.96+np.sqrt(3.32 + 94.6788*wtb -
           (wtb/0.08344)**2)) - 2.7545*wtb*wtb + 5.1799*wtb - 0.2683)

    F_VH = 0.25*np.sin(5.236*wtb) + 0.325

    # Torsion stiffness from assumed torsion frequency, then torque box cross-sectional area, then mass
    GJ = (_OMEGA_T*Omega)**2 * 0.25 * (0.9*b_w) * mtip * r_pylon*r_pylon
    A_tb = 4*GJ / (_G_TB*F_T*t_w*t_w)
    mass_box = A_tb * _DENS_TB * b_w / _E_TB_EFF
    wt_box = mass_box * _GRAV_FTS2  # noqa: F841 (kept for parity with source; not consumed further)

    # Bending stiffnesses from assumed mode frequencies, then spar area, then mass
    EI_C = (_OMEGA_C*Omega)**2 * (1/48.0) * b_w**3 * mtip * f_mode
    EI_B = (_OMEGA_B*Omega)**2 * (1/48.0) * b_w**3 * mtip * f_mode

    EI_Ctb = _E_TB * F_C * A_tb * c_tb*c_tb * 0.25
    EI_Csp = EI_C - EI_Ctb
    EI_Csp = 0 if EI_Csp < 0 else EI_Csp
    A_Csp = EI_Csp / (_E_SP * c_tb*c_tb * 0.25)

    EI_Btb = _E_TB * F_B * A_tb * t_w*t_w * 0.25
    EI_VH  = _E_SP * F_VH * A_Csp * t_w*t_w * 0.25
    EI_Bsp = EI_B - EI_Btb - EI_VH
    EI_Bsp = 0 if EI_Bsp < 0 else EI_Bsp
    A_Bsp = EI_Bsp / (_E_SP * t_w*t_w * 0.25)

    A_sp = A_Csp + A_Bsp
    mass_spar = _CT * A_sp * _DENS_SPAR * b_w / _E_SP_EFF
    mass_prim = mass_box + mass_spar

    wt_spar = mass_spar * _GRAV_FTS2  # noqa: F841 (kept for parity with source; not consumed further)
    wt_prim = mass_prim * _GRAV_FTS2

    # Aerodynamic surfaces: fairings and flaps
    wt_fair = S_fair * _LOADING_FAIR
    wt_flap = S_flap * _LOADING_FLAP

    mass_flap = wt_flap / _GRAV_FTS2
    mass_fair = wt_fair / _GRAV_FTS2

    # Fittings: fraction of everything else so far
    mass_fit = _F_FIT/(1.0-_F_FIT) * (mass_prim+mass_fair+mass_flap)
    wt_fit = mass_fit * _GRAV_FTS2

    # Folding mechanism: fraction of everything else (fixed at 0 — no fold modelled)
    mass_fold = _F_FOLD * (mass_prim+mass_fair+mass_flap+mass_fit)
    wt_fold = mass_fold * _GRAV_FTS2

    # Tilt mechanism: 13% of everything that tilts (tiltwing: whole wing; tiltrotor: tip mass only)
    if tilt_wing:
        wt_tilt = _F_TILT * (wt_prim+wt_fair+wt_flap+wt_fit+wtip)
    else:
        wt_tilt = _F_TILT * wtip

    # Apply technology factor and wing count
    wt_prim = wt_prim * tech_factor * n_wings
    wt_fair = wt_fair * tech_factor * n_wings
    wt_flap = wt_flap * tech_factor * n_wings
    wt_fit  = wt_fit  * tech_factor * n_wings
    wt_fold = wt_fold * tech_factor * n_wings
    wt_tilt = wt_tilt * tech_factor * n_wings

    total = wt_prim + wt_fair + wt_flap + wt_fit + wt_fold + wt_tilt

    return {
        'primary': wt_prim * Units.lb,
        'fairing': wt_fair * Units.lb,
        'flaps':   wt_flap * Units.lb,
        'fitting': wt_fit  * Units.lb,
        'tilt':    wt_tilt * Units.lb,
        'folding': wt_fold * Units.lb,
        'total':   total   * Units.lb,
    }


def compute_wing_tip_mass(wing, vehicle):
    """ Sums the NDARC `Wtip` component group for one wing (NDARC Theory Manual v1.11, pg. 264,
        eq. following 12852): rotor group + engine/nacelle group + drive system + conversion
        flight controls associated with that wing's tilting propulsor group.

        NDARC defines `Wtip` by **propulsor-group membership, not spanwise position** — a motor
        mounted mid-span still belongs to `Wtip` if it tilts with that wing's conversion system.
        The manual's own example of the opposite case ("engine and transmission are not at the tip
        location with the rotor") is handled by an explicit adjustment, not by inferring position
        from geometry. Accordingly this function requires the propulsor-to-wing association to be
        stated explicitly on the wing, not inferred from `motor.origin`:

            wing.NDARC.tip_propulsor_tags   list of propulsor tags whose rotor+motor+hub+servo
                                             mass group with this wing's tip weight. Required.
            wing.NDARC.tip_weight_factor    fWtip, multiplicative [-], default 1.0
            wing.NDARC.tip_weight_increment xWtip, additive [kg], default 0.0

        `fWtip`/`xWtip` match NDARC Input Dictionary v1.19 pg. 143 exactly (`Wtip_adjusted =
        fWtip*Wtip_computed + xWtip`) — the manual gives no formula for either, only a qualitative
        rule ("negative increment required when engine and transmission not at tip location with
        rotor"); sizing the correction is left to the analyst's own engineering judgement of how
        much of the summed component mass is actually displaced away from the tip.

        Only rotor and motor mass are summed for now (both already populated elsewhere in the OEW
        rollup, regardless of which method computed them). Conversion flight-controls mass is not
        yet included — that component (`afdd/flight_controls.py`) is next in the porting order
        (see `20-rcaide-weight-method-porting-inventory.md`) and will be added here once ported,
        not before.

        Inputs:
            wing        RCAIDE Wing Data Structure, with `wing.NDARC.tip_propulsor_tags` set
            vehicle     RCAIDE Vehicle Data Structure

        Outputs:
            tip_mass:   Sum of tagged propulsors' rotor + motor mass, plus adjustment    [kg]
    """
    tip_tags = wing.NDARC.tip_propulsor_tags
    f_wtip = getattr(wing.NDARC, 'tip_weight_factor', 1.0)
    x_wtip = getattr(wing.NDARC, 'tip_weight_increment', 0.0)

    tip_mass = 0.0
    matched_tags = set()
    for network in vehicle.networks:
        for propulsor in network.propulsors:
            if propulsor.tag in tip_tags:
                matched_tags.add(propulsor.tag)
                tip_mass += propulsor.rotor.mass_properties.mass
                tip_mass += propulsor.motor.mass_properties.mass

    missing = set(tip_tags) - matched_tags
    if missing:
        raise ValueError(
            f"wing '{wing.tag}': NDARC.tip_propulsor_tags references propulsor(s) not found on "
            f"the vehicle: {sorted(missing)}"
        )

    return f_wtip * tip_mass + x_wtip
