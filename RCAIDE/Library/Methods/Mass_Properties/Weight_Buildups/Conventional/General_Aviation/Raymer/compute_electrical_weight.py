# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Conventional/General_Aviation/Raymer/compute_electrical_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

from RCAIDE.Framework.Core import Units

# ----------------------------------------------------------------------------------------------------------------------
# Electrical Weight — Raymer (second GA variant, Welec,2)
# ----------------------------------------------------------------------------------------------------------------------
def compute_electrical_weight(vehicle):
    """ Calculates electrical system mass using Raymer's second GA electrical method (Welec,2).

        Source:
            Raymer, Aircraft Design: A Conceptual Approach, as transcribed in Chakraborty &
            Mishra, Journal of Aircraft, Vol. 58, No. 4, 2021, Appendix C, Table C1. Confirmed
            2026-07-28 by direct visual read of the paper's equation image. **Not** independently
            located in the real Raymer textbook text this project has access to (that text's GA
            weights section, eq. 15.47-15.59, only has the `12.57*(Wfs+Wavio)^0.51` variant —
            already implemented as `compute_systems_weight.py`'s `W_electrical` — and the
            Transport-class Raymer electrical formula, `7.291*Rkva^0.782*(2L)^0.346*Neng^0.1`, is a
            structurally different formula entirely, confirmed present in
            `Conventional/Transport/Raymer/compute_systems_weight.py`). This specific
            `0.0268*Wto` variant may come from a different chapter/edition of Raymer's book not
            covered by the available text extraction — treat as sourced from the paper only,
            same confidence level as the Roskam ports, unlike every other Raymer GA formula in
            this folder (all independently confirmed against the primary textbook).

        Welec,2 = 0.0268*Wto

        Inputs:
            vehicle    RCAIDE Vehicle Data Structure

        Outputs:
            weight:    electrical system mass                                             [kg]
    """
    W_to = vehicle.mass_properties.max_takeoff / Units.lb
    weight_english = 0.0268 * W_to
    weight = weight_english * Units.lbs
    return weight
