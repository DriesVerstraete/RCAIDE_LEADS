# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/CADDEE_LPC/compute_empennage_weight.py
#
#

from ._regression import evaluate_empennage_regression

# ----------------------------------------------------------------------------------------------------------------------
#  Constants — nasa_lpc empennage_mass_coeff
# ----------------------------------------------------------------------------------------------------------------------
# Direct port of `m4_nasa_lpc.py`'s `empennage_mass_coeff`. Used live in CADDEE_alpha's own
# reference case (`ex_lpc.py`), with their 1.3x scaler applied. Unlike wing/fuselage/boom, this
# regression is fit against h-tail/v-tail area only (3 coefficients, not the shared 5-variable
# form) — combined h-tail+v-tail mass, matching source (no separate h-tail/v-tail mass split;
# source only splits CG, not mass, 65/35 between them — not ported here, mass only).

_EMPENNAGE_MASS_COEFF = [8.43266623, 10.05410839, -0.19944806479469435]


def compute_empennage_weight(h_tail_area, v_tail_area, correction_factor=1.3):
    """ Calculates combined horizontal+vertical tail mass using CADDEE_alpha's nasa_lpc
        linear-regression surrogate, fit against a DOE sweep of the NASA Lift+Cruise reference
        vehicle's `aframe`-FEA-optimized empennage structure.

        **Valid only for lift+cruise-class vehicles near the LPC reference vehicle's own
        design-space neighborhood** — see module-level docstring in `_regression.py`.

        Source:
            LSDOlab/UCSD `CADDEE_alpha/core/aircraft/models/weights/nasa_lpc/
            m4_nasa_lpc.py::compute_empennage_mps`, `empennage_mass_coeff`.

        Inputs:
            h_tail_area        horizontal tail reference area                          [m^2]
            v_tail_area        vertical tail reference area                            [m^2]
            correction_factor  multiplicative correction, default 1.3 (CADDEE_alpha's own
                                value, applied live in their own reference case for this
                                component; user-overridable)                            [Unitless]

        Outputs:
            mass:   combined horizontal + vertical tail mass                            [kg]
    """
    mass = evaluate_empennage_regression(h_tail_area, v_tail_area, _EMPENNAGE_MASS_COEFF)
    return mass * correction_factor
