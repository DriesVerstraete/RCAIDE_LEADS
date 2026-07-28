# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Conventional/General_Aviation/Roskam/compute_wing_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

import numpy as np
from RCAIDE.Framework.Core import Units

# ----------------------------------------------------------------------------------------------------------------------
# Wing Weight — Roskam (two independent methods)
# ----------------------------------------------------------------------------------------------------------------------
# Both direct ports of Roskam, Airplane Design Parts I Through VIII, 2nd ed., DARcorp., 2003, as
# transcribed in Chakraborty & Mishra, Journal of Aircraft, Vol. 58, No. 4, 2021, Appendix C,
# Table C1 — confirmed 2026-07-28 by direct visual read of the paper's equation image (both text
# extractors mangled these formulas' nested fractions). Not independently cross-checked against
# Roskam's own textbook (not available locally, unlike Raymer's, which was).
#
# `KTW` (tilt-mechanism weight penalty, default 1.20 per the paper) is exposed as a parameter, not
# hardcoded — Chakraborty2022 applies it only for their tiltwing use case; a conventional
# fixed-wing GA aircraft should call with `KTW=1.0`.


def compute_wing_weight_1(vehicle, wing, KTW=1.0):
    """ Calculates wing mass using Roskam's first GA wing method (Wwing,1).

        Wwing,1 = KTW * 96.948 * [ (Wto*nult/1e5)^0.65 * (ARw/cos(Λc/4,w))^0.57 * (Sw/100)^0.61
                                    * ((1+λ)/(2*(t/c)_w))^0.36 * (1+Vmax/500)^0.5 ]^0.993

        Inputs:
            vehicle    RCAIDE Vehicle Data Structure
            wing       RCAIDE Wing Data Structure (main wing)
            KTW        tilt-mechanism weight penalty factor, default 1.0 (Chakraborty2022 uses
                        1.20 for their tiltwing configuration)                          [Unitless]

        Outputs:
            weight:    wing mass                                                          [kg]
    """
    Nult = vehicle.flight_envelope.ultimate_load
    TOW = vehicle.mass_properties.max_takeoff
    Vmax = vehicle.flight_envelope.maneuver.equivalent_speed.velocity_max_dive

    AR_w = wing.aspect_ratio
    taper_w = wing.taper
    t_c_w = wing.thickness_to_chord
    sweep_w = wing.sweeps.quarter_chord

    W_to = TOW / Units.lb
    S_w = wing.areas.reference / Units.ft**2
    V_max = Vmax / Units.kts

    weight_english = KTW * 96.948 * (
        (W_to*Nult/1e5)**.65 * (AR_w/np.cos(sweep_w))**.57 * (S_w/100.)**.61 *
        ((1.+taper_w)/(2.*t_c_w))**.36 * (1.+V_max/500.)**.5
    )**.993

    weight = weight_english * Units.lbs
    return weight


def compute_wing_weight_2(vehicle, wing, sweep_half_chord=None, KTW=1.0):
    """ Calculates wing mass using Roskam's second GA wing method (Wwing,2).

        Wwing,2 = KTW * 0.00125*Wto * (bw/cos(Λc/2,w))^0.75 * (1+6.3*cos(Λc/2,w)/bw)^0.5
                  * (bw*Sw/(tr,w*Wto*cos(Λc/2)))^0.30 * nult^0.55

        Inputs:
            vehicle             RCAIDE Vehicle Data Structure
            wing                RCAIDE Wing Data Structure (main wing)
            sweep_half_chord    wing half-chord sweep — RCAIDE's Wing object does not carry this
                                 natively (only quarter-chord/leading-edge sweep), so it must be
                                 supplied explicitly. Defaults to the wing's quarter-chord sweep if
                                 not given — an approximation, not an exact conversion.       [rad]
            KTW                 tilt-mechanism weight penalty factor, default 1.0
                                 (Chakraborty2022 uses 1.20 for their tiltwing configuration) [Unitless]

        Outputs:
            weight:             wing mass                                                    [kg]
    """
    Nult = vehicle.flight_envelope.ultimate_load
    TOW = vehicle.mass_properties.max_takeoff

    if sweep_half_chord is None:
        sweep_half_chord = wing.sweeps.quarter_chord

    b_w = wing.spans.projected / Units.ft
    S_w = wing.areas.reference / Units.ft**2
    t_r_w = (wing.thickness_to_chord * wing.chords.root) / Units.ft
    W_to = TOW / Units.lb

    weight_english = KTW * 0.00125 * W_to * (b_w/np.cos(sweep_half_chord))**.75 * \
        (1.+6.3*np.cos(sweep_half_chord)/b_w)**.5 * \
        ((b_w*S_w)/(t_r_w*W_to*np.cos(sweep_half_chord)))**.30 * Nult**.55

    weight = weight_english * Units.lbs
    return weight
