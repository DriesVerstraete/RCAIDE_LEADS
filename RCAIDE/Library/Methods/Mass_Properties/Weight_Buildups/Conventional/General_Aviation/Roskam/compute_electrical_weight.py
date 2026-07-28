# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Conventional/General_Aviation/Roskam/compute_electrical_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

from RCAIDE.Framework.Core import Units

# ----------------------------------------------------------------------------------------------------------------------
# Electrical Weight — Roskam
# ----------------------------------------------------------------------------------------------------------------------
def compute_electrical_weight(fuel_system_weight, avionics_weight):
    """ Calculates electrical system mass using Roskam's GA method (Welec,1).

        Source:
            Roskam, Airplane Design Parts I Through VIII, 2nd ed., DARcorp., 2003 — as
            transcribed in Chakraborty & Mishra, Journal of Aircraft, Vol. 58, No. 4, 2021,
            Appendix C, Table C1. Confirmed 2026-07-28 by direct visual read of the paper's
            equation image.

            Note: this exact formula (coefficient 12.57, exponent 0.51) is ALSO confirmed present
            verbatim in the real Raymer textbook (eq. 15.53, checked directly 2026-07-28 against
            the parsed textbook text) — RCAIDE's existing `Raymer/compute_systems_weight.py`
            already implements it, labeled Raymer, not Roskam. Both citations may be genuine (the
            two textbooks may share the regression) — not resolved either way, since only
            Raymer's textbook was available locally to check. Implemented as its own function
            here regardless, matching the paper's own Table C1 attribution, so it can be called
            independently of Raymer's bundled `compute_systems_weight`.

        Welec,1 = 12.57*(Wfs+Wavio)^0.51

        Inputs:
            fuel_system_weight    fuel system mass                                        [kg]
            avionics_weight       avionics mass                                            [kg]

        Outputs:
            weight:                electrical system mass                                  [kg]
    """
    W_fs = fuel_system_weight / Units.lb
    W_avio = avionics_weight / Units.lb

    weight_english = 12.57 * (W_fs + W_avio)**.51

    weight = weight_english * Units.lbs
    return weight
