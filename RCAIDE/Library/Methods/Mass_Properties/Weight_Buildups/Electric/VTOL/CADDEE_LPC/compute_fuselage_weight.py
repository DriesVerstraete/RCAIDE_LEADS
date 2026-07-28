# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/CADDEE_LPC/compute_fuselage_weight.py
#
#

from ._regression import evaluate_regression

# ----------------------------------------------------------------------------------------------------------------------
#  Constants — nasa_lpc fuselage_reg coefficients
# ----------------------------------------------------------------------------------------------------------------------
# Direct port of `m4_nasa_lpc.py`'s `fuselage_reg['fuselage_mass']` coefficients. Unlike wing, this
# IS the regression CADDEE_alpha's own reference case (`ex_lpc.py`) actually uses live, always with
# their own 1.3x `scaler` applied (their comment: "M4-regression mass models (scaled to match
# better with NDARC)"). Cross-checked against the real `lift_cruise` vehicle and real NDARC AFDD84
# fuselage weight on the same vehicle, 2026-07-28: raw regression = 62% of NDARC, x1.3 = 81% of
# NDARC, Vahana = 48% of NDARC. The 1.3x claim is directionally true (closer to NDARC) but not an
# exact match on this vehicle — consistent with the scaler having been hand-tuned against NASA's
# own LPC vehicle specifically. See
# 01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md for the full
# comparison table.

_FUSELAGE_MASS_COEFF = [1.01472161e+00, -4.06758251e-01, 4.25974124e+01, 3.10575276e-02,
                         6.87345416e-02, -1.16727769e+02]


def compute_fuselage_weight(wing_area, wing_AR, fuselage_length, battery_mass, cruise_speed,
                             correction_factor=1.3):
    """ Calculates fuselage mass using CADDEE_alpha's nasa_lpc linear-regression surrogate, fit
        against a DOE sweep of the NASA Lift+Cruise reference vehicle's `aframe`-FEA-optimized
        fuselage structure. This is the one nasa_lpc regression CADDEE_alpha's own reference case
        actually uses live (with their own 1.3x scaler) — see module docstring.

        **Valid only for lift+cruise-class vehicles near the LPC reference vehicle's own
        design-space neighborhood** — see module-level docstring in `_regression.py`.

        Source:
            LSDOlab/UCSD `CADDEE_alpha/core/aircraft/models/weights/nasa_lpc/
            m4_nasa_lpc.py::compute_fuselage_mps`, `fuselage_reg['fuselage_mass']` coefficients.

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
            mass:   fuselage mass                                                       [kg]
    """
    mass = evaluate_regression(wing_area, wing_AR, fuselage_length, battery_mass, cruise_speed,
                                _FUSELAGE_MASS_COEFF)
    return mass * correction_factor
