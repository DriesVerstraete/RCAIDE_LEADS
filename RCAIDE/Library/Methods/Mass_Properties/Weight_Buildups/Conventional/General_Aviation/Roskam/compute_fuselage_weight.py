# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Conventional/General_Aviation/Roskam/compute_fuselage_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

from RCAIDE.Framework.Core import Units

# ----------------------------------------------------------------------------------------------------------------------
# Fuselage Weight — Roskam
# ----------------------------------------------------------------------------------------------------------------------
def compute_fuselage_weight(fuselage, vehicle):
    """ Calculate the weight of a fuselage for a GA aircraft using Roskam's method.

        Source:
            Roskam, J., Airplane Design Parts I Through VIII, 2nd ed., DARcorp., 2003 — as
            transcribed in Chakraborty, I., and Mishra, A. A., "Generalized Energy-Based Flight
            Vehicle Sizing and Performance Analysis Methodology," Journal of Aircraft, Vol. 58,
            No. 4, 2021, Appendix C, Table C1, Eq. Wfus,1. Confirmed 2026-07-28 by direct visual
            read of the paper's equation image (both text-extraction tools mangled this formula's
            nested-fraction layout) — not independently cross-checked against Roskam's own
            textbook (not available locally, unlike Raymer's, which was).

        Wfus,1 = 200 * [ (Wto*nult/1e5)^0.286 * (Lfus/10)^0.857 * ((wfus+hfus)/10) * (Vc/100)^0.338 ]^1.1

        Inputs:
            fuselage    RCAIDE Fuselage Data Structure
            vehicle     RCAIDE Vehicle Data Structure

        Outputs:
            weight:     fuselage mass                                                    [kg]
    """
    Nult = vehicle.flight_envelope.ultimate_load
    TOW = vehicle.mass_properties.max_takeoff
    Vc = vehicle.flight_envelope.maneuver.equivalent_speed.velocity_max_cruise

    L_fus = fuselage.lengths.total / Units.ft
    w_fus = fuselage.width / Units.ft
    h_fus = fuselage.heights.maximum / Units.ft
    W_to = TOW / Units.lb
    V_c = Vc / Units.kts

    weight_english = 200. * ((W_to*Nult/1e5)**.286 * (L_fus/10.)**.857 *
                              ((w_fus+h_fus)/10.) * (V_c/100.)**.338)**1.1

    weight = weight_english * Units.lbs
    return weight
