# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/NDARC/compute_rotor_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

# RCAIDE imports
from RCAIDE.Framework.Core import Units

# ----------------------------------------------------------------------------------------------------------------------
#  Constants
# ----------------------------------------------------------------------------------------------------------------------
# Direct port of Hydra's `src/Python/Stage_1/afdd/rotor_wt.py::rotor_weight`
# (`github.com/VahanaOpenSource/vtol_sizing`) — blades, hub/hinge, spinner, and blade-fold
# structure. Source is entirely imperial-calibrated (its own header: "ALL UNITS IN IMPERIAL
# (FPS)"); SI inputs converted internally via RCAIDE's own `Units` module. See
# 01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md.

_M2F = 1.0 / Units.ft
_KG2LB = 1.0 / Units.lb

_F_FOLD = 0.0   # blade fold weight fraction — no fold mechanism modelled (fixed at 0 in source too)


def compute_rotor_weight(aircraft_id, n_blade, n_rotor, radius, chord, tip_speed, nu_blade,
                          tech_factor=1.0):
    """ Calculates rotor blade, hub, spinner, and blade-fold structure mass using either the
        AFDD82 or AFDD00 statistical regression, selected by `aircraft_id`.

        Model selection (matches source exactly, including a quirk found during translation):
        `aircraft_id == 2` (tiltrotor) selects the AFDD00 model, with the tilt-factored blade
        weight (`f_tilt=1.1794`) and a spinner sized from `dia_spin = 0.2*R`. Every other
        `aircraft_id` value (1, 3, 5, 6, 7) selects the AFDD82 model, `f_tilt=1.0`, no spinner.

        **Known source quirk, ported faithfully rather than "fixed"**: the AFDD00 block contains a
        coaxial-specific hub formula gated on `aircraft_id == 3`, citing Johnson/Moodie/Yeo
        "Design and Performance of Lift-Offset Rotorcraft for Short-Haul Missions" — but the outer
        model-selection logic routes `aircraft_id == 3` to AFDD82, not AFDD00, so that coaxial hub
        formula is unreachable in the original source. Replicated here as unreachable too (an
        `aircraft_id == 3` case still runs AFDD82, exactly like the source), not "corrected" to
        make the coaxial formula reachable — a real port should reproduce the source's actual
        behavior, not its apparent intent. Flagged in case reachability is wanted later.

        The source's unrecognized-`aircraft_id` error path calls `sys.exit(1)` without importing
        `sys` (would raise `NameError`, not exit cleanly, if ever triggered) — replaced here with a
        proper `ValueError`, a deliberate deviation from source (bug in source's error handling,
        not core logic).

        Source:
            Hydra `afdd/rotor_wt.py::rotor_weight`, AFDD82/AFDD00 models,
            NDARC Theory Manual v1.11 Section 29 (rotor group).

        Inputs:
            aircraft_id     Hydra vehicle-classification flag — 2 selects AFDD00 (tiltrotor
                             tilt-factored), all other recognized values (1,3,5,6,7) select AFDD82  [int]
            n_blade         number of blades per rotor                                      [Unitless]
            n_rotor         number of rotors in this group                                  [Unitless]
            radius          rotor radius                                                    [m]
            chord           blade chord                                                     [m]
            tip_speed       rotor tip speed                                                 [m/s]
            nu_blade        blade flap natural frequency ratio (per-rev)                    [Unitless]
            tech_factor     technology weight-scaling factor                                [Unitless]

        Outputs:
            weight:            dict with 'blades', 'hub', 'spinner', 'folding', all           [kg]
            mass_kg_per_unit:  average per-rotor mass (blades+hub+spinner+folding)/n_rotor     [kg]
    """
    R = radius * _M2F
    c = chord * _M2F
    v_tip = tip_speed * _M2F
    nu_hub = nu_blade

    if aircraft_id in (1, 3, 5, 6, 7):
        f_tilt = 1.0
        dia_spin = 0.0
        model = 'afdd82'
    elif aircraft_id == 2:
        f_tilt = 1.1794
        dia_spin = 0.2 * R
        model = 'afdd00'
    else:
        raise ValueError(f"compute_rotor_weight: unrecognized aircraft_id {aircraft_id!r}")

    if model == 'afdd00':
        wght_blade = (0.0024419 * f_tilt *
                      n_rotor * n_blade**0.53479 *
                      R**1.74231 * c**0.77291 *
                      v_tip**0.87562 * nu_blade**2.51048)

        if aircraft_id == 3:
            # Unreachable — see docstring. aircraft_id==3 always routes to AFDD82 above.
            wght_hub = (0.0061182 * n_rotor * n_blade**0.20373 *
                        R**0.60406 * v_tip**0.52803 *
                        nu_hub**1.00218 * (wght_blade/n_rotor)**0.87127)
        else:
            wght_hub = (0.1837 * n_rotor * n_blade**0.16383 *
                        R**0.19937 * v_tip**0.06171 *
                        nu_hub**0.46203 * (wght_blade/n_rotor)**1.02958)
    else:
        wght_blade = (0.02606*n_rotor * n_blade**0.6592 *
                      R**1.3371 * c**0.9959 *
                      v_tip**0.6682 * nu_blade**2.5279)

        wght_hub = (0.003722*n_rotor * n_blade**0.2807 *
                    R**1.5377 * v_tip**0.4290 *
                    nu_hub**2.1414 * (wght_blade/n_rotor)**0.5505)

    wght_spin = 7.386*n_rotor*dia_spin**2
    wght_fold = _F_FOLD*wght_blade

    wght_blade = wght_blade * tech_factor
    wght_hub = wght_hub * tech_factor
    wght_spin = wght_spin * tech_factor
    wght_fold = wght_fold * tech_factor

    weight = {
        'blades': wght_blade * Units.lb,
        'hub': wght_hub * Units.lb,
        'spinner': wght_spin * Units.lb,
        'folding': wght_fold * Units.lb,
    }

    mass_kg_per_unit = (wght_blade + wght_hub + wght_spin + wght_fold) / n_rotor * Units.lb

    return weight, mass_kg_per_unit
