# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/NDARC/compute_flight_control_system_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

# RCAIDE imports
from RCAIDE.Framework.Core import Units

# ----------------------------------------------------------------------------------------------------------------------
#  Constants — AFDD82 flight control system model
# ----------------------------------------------------------------------------------------------------------------------
# Direct port of Hydra's `src/Python/Stage_1/afdd/flight_controls.py::flightctrl_weight`
# (`github.com/VahanaOpenSource/vtol_sizing`). Covers rotary-wing flight controls, fixed-wing
# flight controls, conversion (tilt) boost mechanism, and hydraulics for all of the above — the
# source does not separate "tilting hardware" out as its own function; the conversion-mechanism
# terms (`wght_CVmb`/`wght_CVnb`/`wght_CVhyd`, NDARC's biggest structural gap in RCAIDE) are
# threaded through this file's combined rotary-wing/hydraulics totals via `aircraftID==2`
# ("tilt-rotor" per `afdd/rotor_wt.py`'s own comment), not isolable without restructuring the
# source's own combination logic. Ported whole for fidelity; the individual terms are still
# returned separately in the output dict.
#
# All inputs are in the imperial units this AFDD82 regression was calibrated in (ft, lb, ft/s);
# SI inputs are converted internally via RCAIDE's own `Units` module. See
# 01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md.

_M2F = 1.0 / Units.ft
_KG2LB = 1.0 / Units.lb

_F_RWNB  = 1.25     # fraction rotary-wing non-boosted weight
_F_RWHYD = 0.4      # fraction rotary-wing hydraulic weight
_F_RWRED = 3.0      # redundancy factor
_F_CVNB  = 0.1      # fraction conversion non-boosted weight (of mass_CVmb)
_F_CVMB  = 0.02     # fraction conversion boosted mechanism weight (of MTOW)
_F_MBSV  = 1.3029   # ballistic survivability factor, boosted
_F_BSV   = 1.117    # ballistic survivability factor, non-boosted
_F_CVHYD = 0.1      # fraction conversion hydraulic weight (of boost mechanism)


def compute_flight_control_system_weight(n_rotor, n_blade, chord, tip_speed, vehicle_mtow,
                                           aircraft_id, n_prop, n_wing, rotor_radius,
                                           tech_factor=1.0):
    """ Calculates the mass of the flight control system: rotary-wing controls, fixed-wing
        controls, conversion (tilt) boost mechanism, and hydraulics for all of the above — the
        AFDD82 statistical regression NDARC uses (Theory Manual v1.11 Section 29-8/29-9).

        This is NDARC's biggest structural gap versus RCAIDE's native `Physics_Based` method,
        which does not model flight control system or tilting-hardware weight at all
        (`output.empty.systems.control_systems` is hardcoded 0.0).

        The conversion boost mechanism (`aircraft_id == 2`, i.e. tiltrotor per the source's own
        convention in `afdd/rotor_wt.py`) is the tilting-hardware term specifically:
        `wght_CVmb = f_CVmb * MTOW`, `wght_CVnb = f_CVnb * wght_CVmb`, `wght_CVhyd = f_CVhyd *
        wght_CVmb`. It is NOT active for any other `aircraft_id` value.

        Assumptions:
            `aircraft_id` follows Hydra's own internal vehicle-classification convention, not
            RCAIDE's. Confirmed meanings: 2 = tiltrotor (activates conversion controls). 5 =
            disables the fixed-wing actuator weight term entirely (meaning unconfirmed — not yet
            traced to a specific vehicle archetype in the source). All other values: standard
            rotary + fixed-wing control weight, no conversion mechanism. Mapping a RCAIDE vehicle's
            architecture onto this convention (e.g. "does this vehicle have tilting propulsors")
            is deliberately left to the caller, not inferred here — same reasoning as
            `compute_wing_tip_mass`'s explicit-tagging approach.

        Source:
            Hydra `afdd/flight_controls.py::flightctrl_weight`, AFDD82 model,
            NDARC Theory Manual v1.11 Section 29-8/29-9.

        Inputs:
            n_rotor         number of primary (main) rotors                      [Unitless]
            n_blade         number of blades per primary rotor                   [Unitless]
            chord           primary rotor blade chord                           [m]
            tip_speed       primary rotor tip speed                              [m/s]
            vehicle_mtow    vehicle max takeoff weight                           [kg]
            aircraft_id     Hydra vehicle-classification flag (see Assumptions)  [int]
            n_prop          number of auxiliary thrust propellers (0 if none)    [Unitless]
            n_wing          number of wings (0 if none — tailless/no-wing case)  [Unitless]
            rotor_radius    primary rotor radius (used for the no-wing Htail-only
                             actuator sizing branch; only consumed when n_wing==0) [m]
            tech_factor     technology weight-scaling factor                     [Unitless]

        Outputs:
            weight:   dict with 'hydraulics', 'rw_flt_ctrl', 'fw_flt_ctrl', 'total', all   [kg]
    """
    chord_ft = chord * _M2F
    v_tip = tip_speed * _M2F
    wmto = vehicle_mtow * _KG2LB
    r_ft = rotor_radius * _M2F

    if n_wing > 0:
        w = 0.91 * (wmto**0.6)
    else:
        s_ht = 45.0 * (r_ft / 26.83)**2
        w = 0.01735 * (wmto**0.64345) * (s_ht**0.40952)

    w_fc = (0.2873 * _F_MBSV * (n_rotor * n_blade)**0.6257 *
            chord_ft**1.3286 * (0.01*v_tip)**2.1129 *
            _F_RWRED**0.8942)

    wght_RWb = (0.02324 * _F_BSV * (n_rotor*n_blade)**1.0042 *
                n_rotor**0.1155 * chord_ft**2.2296 *
                (0.01*v_tip)**3.1877)

    if n_prop > 0:
        v_tip_prop = 600.0
        chord_prop = chord_ft*0.2
        n_rotor_prop = n_prop
        n_blade_prop = 4

        w_fc = w_fc + (0.2873 * _F_MBSV * (n_rotor_prop * n_blade_prop)**0.6257 *
                       chord_prop**1.3286 * (0.01*v_tip_prop)**2.1129 *
                       _F_RWRED**0.8942)

        wght_RWb = wght_RWb + (0.02324 * _F_BSV * (n_rotor_prop*n_blade_prop)**1.0042 *
                                n_rotor_prop**0.1155 * chord_prop**2.2296 *
                                (0.01*v_tip_prop)**3.1877)

    wght_RWmb = (1-_F_RWHYD)*w_fc
    wght_RWnb = _F_RWNB*(1-_F_RWHYD)*w_fc

    if aircraft_id == 2:
        wght_CVmb = _F_CVMB*wmto
        wght_CVnb = _F_CVNB*wght_CVmb
    else:
        wght_CVmb = 0.0
        wght_CVnb = 0.0

    wght_avionics = 0.0

    all_fltcon = wght_RWmb + wght_RWb + wght_RWnb + wght_CVmb + wght_CVnb + wght_avionics

    if aircraft_id == 5:
        fx_wing_ctrl = 0.0
    else:
        fx_wing_ctrl = w

    rt_wing_ctrl = all_fltcon

    wght_RWhyd = _F_RWHYD*w_fc
    wght_CVhyd = _F_CVHYD*wght_CVmb

    wght_hydraulics = wght_RWhyd + wght_CVhyd

    wght_hydraulics = wght_hydraulics * tech_factor
    rt_wing_ctrl = rt_wing_ctrl * tech_factor
    fx_wing_ctrl = fx_wing_ctrl * tech_factor

    total = wght_hydraulics + rt_wing_ctrl + fx_wing_ctrl

    return {
        'hydraulics': wght_hydraulics * Units.lb,
        'rw_flt_ctrl': all_fltcon * Units.lb,
        'fw_flt_ctrl': fx_wing_ctrl * Units.lb,
        'total': total * Units.lb,
    }
