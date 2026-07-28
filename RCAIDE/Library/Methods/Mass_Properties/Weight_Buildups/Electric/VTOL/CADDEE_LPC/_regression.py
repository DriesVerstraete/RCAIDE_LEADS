# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/CADDEE_LPC/_regression.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  Shared regression evaluator — nasa_lpc linear-regression coefficient sets
# ----------------------------------------------------------------------------------------------------------------------
# Direct port of LSDOlab/UCSD's `CADDEE_alpha/core/aircraft/models/weights/nasa_lpc/
# m4_nasa_lpc.py::evaluate_regression`/`evaluate_empennage_regression` (`github.com/LSDOlab/
# CADDEE_alpha`). Coefficients are a linear-regression fit against outputs of a full FEA
# structural model (`aframe`, MIT-licensed 1D beam FEA) that was itself optimized, run across a
# DOE sweep of the NASA Lift+Cruise (LPC) reference vehicle geometry (Ruh & Gandhi, "A Parametric
# Mass Estimation Method for Electric Vertical Takeoff and Landing Aircraft," AIAA 2023-4280).
#
# **Valid only for lift+cruise-class vehicles, and only near the LPC reference vehicle's own
# design-space neighborhood** — this is a surrogate for "what the optimized FEA structure weighs
# at this point in the design space", not a broad multi-vehicle statistical fit like NDARC's
# AFDD82/84/00 regressions. Do not use for tiltrotor, tiltwing, or any vehicle whose geometry sits
# far outside a conventional L+C layout (fixed wing, cruise pusher, boom-mounted lift rotors).

def evaluate_regression(wing_area, wing_AR, fuselage_length, battery_mass, cruise_speed, coeffs):
    """ Direct port of `m4_nasa_lpc.py::evaluate_regression` — a 6-coefficient linear fit against
        5 design variables, shared by the wing/fuselage/boom mass (and CG/inertia, not ported
        here) regressions.

        qty = c0*wing_area + c1*wing_AR + c2*fuselage_length + c3*battery_mass + c4*cruise_speed + c5
    """
    return (coeffs[0]*wing_area + coeffs[1]*wing_AR + coeffs[2]*fuselage_length +
            coeffs[3]*battery_mass + coeffs[4]*cruise_speed + coeffs[5])


def evaluate_empennage_regression(h_tail_area, v_tail_area, coeffs):
    """ Direct port of `m4_nasa_lpc.py::evaluate_empennage_regression` — a 3-coefficient linear fit
        against h-tail and v-tail area only (empennage mass is fit separately from the other
        components, not part of the 5-variable regression above).

        qty = c0*h_tail_area + c1*v_tail_area + c2
    """
    return coeffs[0]*h_tail_area + coeffs[1]*v_tail_area + coeffs[2]
