# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Conventional/General_Aviation/Roskam/compute_vertical_tail_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

from RCAIDE.Framework.Core import Units

# ----------------------------------------------------------------------------------------------------------------------
# Vertical Tail Weight — Roskam
# ----------------------------------------------------------------------------------------------------------------------
def compute_vertical_tail_weight(wing, vehicle, n_vertical_tails=1):
    """ Calculates vertical tail mass using Roskam's GA method (Wvt,1).

        Source:
            Roskam, Airplane Design Parts I Through VIII, 2nd ed., DARcorp., 2003 — as
            transcribed in Chakraborty & Mishra, Journal of Aircraft, Vol. 58, No. 4, 2021,
            Appendix C, Table C1. Confirmed 2026-07-28 by direct visual read of the paper's
            equation image. Not independently cross-checked against Roskam's own textbook (not
            available locally, unlike Raymer's, which was).

        Wvt,1 = 98.5*Nvt*[0.289*(Wto*nult/1e5)^0.87*(Svt/100)^1.2*(bvt/λvt)^0.5]^0.458

        Inputs:
            wing                RCAIDE Wing Data Structure (vertical tail)
            vehicle             RCAIDE Vehicle Data Structure
            n_vertical_tails    number of vertical tails, default 1                    [Unitless]

        Outputs:
            weight:             vertical tail mass                                       [kg]
    """
    Nult = vehicle.flight_envelope.ultimate_load
    TOW = vehicle.mass_properties.max_takeoff

    S_vt = wing.areas.reference / Units.ft**2
    b_vt = wing.spans.projected / Units.ft
    taper_vt = wing.taper
    W_to = TOW / Units.lb

    weight_english = 98.5 * n_vertical_tails * (
        0.289 * (W_to*Nult/1e5)**.87 * (S_vt/100.)**1.2 * (b_vt/taper_vt)**.5
    )**.458

    weight = weight_english * Units.lbs
    return weight
