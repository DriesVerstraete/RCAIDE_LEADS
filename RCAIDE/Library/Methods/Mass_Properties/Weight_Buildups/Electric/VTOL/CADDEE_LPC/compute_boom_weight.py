# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/CADDEE_LPC/compute_boom_weight.py
#
#

from ._regression import evaluate_regression

# ----------------------------------------------------------------------------------------------------------------------
#  Constants — nasa_lpc boom_mass_coeffs
# ----------------------------------------------------------------------------------------------------------------------
# Direct port of `m4_nasa_lpc.py`'s `boom_mass_coeffs`. Used live in CADDEE_alpha's own reference
# case (`ex_lpc.py`), with their 1.3x scaler applied — same "matches better with NDARC" claim as
# fuselage. Total mass covers all 4 booms (left/right x inner/outer) combined; source splits it
# evenly (`total_boom_mass / 4`) for CG purposes, not ported here (mass only).

_BOOM_MASS_COEFF = [3.33883108e+00, 5.40016118e+00, -9.17180256e-02, 1.28646937e-04,
                     -6.25697773e-03, -1.17008833e+01]


def compute_boom_weight(wing_area, wing_AR, fuselage_length, battery_mass, cruise_speed,
                         correction_factor=1.3):
    """ Calculates total boom mass (all 4 booms combined) using CADDEE_alpha's nasa_lpc
        linear-regression surrogate, fit against a DOE sweep of the NASA Lift+Cruise reference
        vehicle's `aframe`-FEA-optimized boom structures.

        **Valid only for lift+cruise-class vehicles near the LPC reference vehicle's own
        design-space neighborhood** — see module-level docstring in `_regression.py`.

        Source:
            LSDOlab/UCSD `CADDEE_alpha/core/aircraft/models/weights/nasa_lpc/
            m4_nasa_lpc.py::compute_boom_mps`, `boom_mass_coeffs`.

        Inputs:
            wing_area          wing reference area                                    [m^2]
            wing_AR            wing aspect ratio                                      [Unitless]
            fuselage_length    fuselage length                                        [m]
            battery_mass       battery mass                                           [kg]
            cruise_speed       cruise speed                                           [m/s]
            correction_factor  multiplicative correction, default 1.3 (CADDEE_alpha's own
                                value, applied live in their own reference case for this
                                component; user-overridable)                            [Unitless]

        Outputs:
            mass:   total boom mass, all 4 booms combined                               [kg]
    """
    mass = evaluate_regression(wing_area, wing_AR, fuselage_length, battery_mass, cruise_speed,
                                _BOOM_MASS_COEFF)
    return mass * correction_factor
