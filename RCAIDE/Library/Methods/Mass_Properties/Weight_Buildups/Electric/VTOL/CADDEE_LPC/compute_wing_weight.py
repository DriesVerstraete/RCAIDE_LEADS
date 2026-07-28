# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/CADDEE_LPC/compute_wing_weight.py
#
#

from ._regression import evaluate_regression

# ----------------------------------------------------------------------------------------------------------------------
#  Constants — nasa_lpc wing_reg coefficients
# ----------------------------------------------------------------------------------------------------------------------
# Direct port of `m4_nasa_lpc.py`'s `wing_reg['wing_mass']` coefficients. CADDEE_alpha's own
# reference case (`examples/advanced_examples/ex_lpc.py`) never actually exercises this regression
# for wing (it uses real `aframe` FEA instead — the `compute_wing_mps` call is present in their
# source but commented out) — so there is no equivalent "does the correction factor work" evidence
# for wing the way there is for fuselage/boom/empennage (see module docstring in
# `compute_fuselage_weight.py`). Cross-checked against the real `lift_cruise` vehicle 2026-07-28:
# raw regression runs 42% heavier than Vahana/NDARC/Hydra's wing weight on that vehicle, 84%
# heavier with the 1.3x correction factor — see
# 01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md for the full
# comparison table (including NDARC's own AFDD93 alternative, which sits closer to a 4-method
# convergence cluster than this regression does).

_WING_MASS_COEFF = [1.11379136e+01, 3.14761829e+01, 7.89132288e-01, -2.14257921e-02,
                     2.40041303e-01, -3.20236992e+02]


def compute_wing_weight(wing_area, wing_AR, fuselage_length, battery_mass, cruise_speed,
                         correction_factor=1.3):
    """ Calculates wing mass using CADDEE_alpha's nasa_lpc linear-regression surrogate, fit
        against a DOE sweep of the NASA Lift+Cruise reference vehicle's `aframe`-FEA-optimized
        wing structure.

        **Valid only for lift+cruise-class vehicles near the LPC reference vehicle's own
        design-space neighborhood** — see module-level docstring in `_regression.py`.

        Source:
            LSDOlab/UCSD `CADDEE_alpha/core/aircraft/models/weights/nasa_lpc/
            m4_nasa_lpc.py::compute_wing_mps`, `wing_reg['wing_mass']` coefficients.

        Inputs:
            wing_area          wing reference area                                    [m^2]
            wing_AR            wing aspect ratio                                      [Unitless]
            fuselage_length    fuselage length                                        [m]
            battery_mass       battery mass                                           [kg]
            cruise_speed       cruise speed                                           [m/s]
            correction_factor  multiplicative correction, default 1.3 (CADDEE_alpha's own
                                "matches better with NDARC" factor — validated by them for
                                fuselage/boom/empennage only, NOT wing; user-overridable)  [Unitless]

        Outputs:
            mass:   wing mass                                                          [kg]
    """
    mass = evaluate_regression(wing_area, wing_AR, fuselage_length, battery_mass, cruise_speed,
                                _WING_MASS_COEFF)
    return mass * correction_factor
