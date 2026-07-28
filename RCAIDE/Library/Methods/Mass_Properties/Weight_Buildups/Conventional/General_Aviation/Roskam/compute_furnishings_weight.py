# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Conventional/General_Aviation/Roskam/compute_furnishings_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

from RCAIDE.Framework.Core import Units

# ----------------------------------------------------------------------------------------------------------------------
# Furnishings Weight — Roskam (two independent methods)
# ----------------------------------------------------------------------------------------------------------------------
# Both direct ports of Roskam, as transcribed in Chakraborty & Mishra, Journal of Aircraft,
# Vol. 58, No. 4, 2021, Appendix C, Table C1. Confirmed 2026-07-28 by direct visual read of the
# paper's equation image. `compute_furnishings_weight_1` is ALSO confirmed present verbatim in
# the real Raymer textbook (eq. 15.59) — RCAIDE's existing `Raymer/compute_systems_weight.py`
# already implements it (labeled Raymer there). Implemented here too, matching the paper's own
# Table C1 attribution, so both Roskam variants can be averaged independently of Raymer's file.


def compute_furnishings_weight_1(vehicle):
    """ Wfur,1 = 0.0582*Wdg - 65

        Inputs:
            vehicle    RCAIDE Vehicle Data Structure

        Outputs:
            weight:    furnishings mass                                                   [kg]
    """
    W_dg = vehicle.mass_properties.max_takeoff / Units.lb
    weight_english = 0.0582*W_dg - 65.
    weight = weight_english * Units.lbs
    return weight


def compute_furnishings_weight_2(vehicle):
    """ Wfur,2 = 0.412*Npax^1.145*Wto^0.489

        Inputs:
            vehicle    RCAIDE Vehicle Data Structure (`number_of_passengers` includes crew,
                        matching the paper's `Npax` definition — "Number of passengers,
                        including crew")

        Outputs:
            weight:    furnishings mass                                                   [kg]
    """
    N_pax = vehicle.number_of_passengers
    W_to = vehicle.mass_properties.max_takeoff / Units.lb
    weight_english = 0.412 * (N_pax**1.145) * (W_to**.489)
    weight = weight_english * Units.lbs
    return weight
