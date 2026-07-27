# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/NDARC/compute_fuselage_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

# RCAIDE imports
from RCAIDE.Framework.Core import Units

# package imports
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  Constants — AFDD84 universal fuselage weight model
# ----------------------------------------------------------------------------------------------------------------------
# Direct port of Hydra's `src/Python/Stage_1/afdd/fuselage.py::fuselage_weight`
# (`github.com/VahanaOpenSource/vtol_sizing`). Source is entirely imperial-calibrated (its own
# header: "ALL UNITS IN IMPERIAL"); SI inputs converted internally via RCAIDE's own `Units` module.
# See 01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md.
#
# Source has two module-level constants (`a=-2.3979`, `b=1.0`) that are defined but never actually
# used in the formula (only `c`/`d` feed the wetted-area regression) — dropped here, not ported,
# since porting an unused constant would be noise. Source also defines `nz` twice at module level
# (4.0, then overwritten to 3.5) — the second value is what's actually used; ported as 3.5 only.
#
# Spot-checked against the real NDARC v1.19 Fortran source (2026-07-27,
# weight_model.f90::GetWeightAFDD_Fuselage, local path:
# `/Users/dverstraete/Sydney Uni Dropbox/Dries Verstraete/Software/NDARC/1_19/NDARC_v1_19_source/`).
# Core AFDD84 regression coefficients (25.41, 0.4879, 0.2075, 0.1676, 0.1512) match exactly — this
# is a faithful port of the main formula.
#
# **Correction, 2026-07-27**: the landing-gear-location/retraction premium (`fLGloc`/`fLGret`) was
# initially hardcoded off here on the (unchecked, wrong) assumption that all current vehicles have
# wing-mounted gear. Checked directly: `lift_cruise/vehicle.py` defines `main_gear`/`nose_gear`
# with fuselage-centerline origins (`[[4.0,0,0]]`/`[[0.5,0,0]]`) — genuinely fuselage-mounted, not
# wing-mounted. `tilt_stopped_rotor_v_tail` and `tiltrotor` define no landing gear at all (zero LG
# mass in any method currently, a separate gap). Now a real switchable input
# (`landing_gear_on_fuselage`), matching real NDARC's `place_LG` conditional
# (`fLGloc=1.1627`/`fLGret=1.1437` when gear is on-fuselage, per `weight_model.f90` lines 329-330).
# `has_cargo_ramp` added for the same reason (`framp`), though no current vehicle has one (still
# defaults off, kept switchable rather than hardcoded for consistency with the LG fix).
#
# One remaining real discrepancy, still dormant, not yet fixed: fold-weight base quantity does NOT
# match NDARC. This port (and Hydra) computes tail-fold weight as a fraction of *fuselage* weight
# and wing-fold weight as a fraction of *(fuselage+tail-fold)* weight. Real NDARC computes
# tail-fold as a fraction of TAIL weight and wing-fold as a fraction of WING weight — a different
# physical basis. Currently dormant (`_F_TFOLD`/`_F_WFOLD` are both 0.0), so numerically inert —
# but wrong if fold weight is ever activated using these fractions. Needs a real fix (passing
# tail/wing weight in) before fold weight can be trusted, not just a coefficient tweak.

_M2F = 1.0 / Units.ft
_KG2LB = 1.0 / Units.lb

_F_LGLOC_ON_FUSELAGE = 1.1627   # landing-gear-on-fuselage location factor (real NDARC value)
_F_LGRET_RETRACTABLE = 1.1437   # retractable landing gear factor (real NDARC value, only applies
                                 # in combination with on-fuselage gear, matching NDARC's own
                                 # `(place_LG<=1) AND (kind_LG>=1)` condition)
_F_RAMP_CARGO        = 1.2749   # cargo ramp factor (real NDARC value)
_F_TFOLD = 0.0    # tail fold weight fraction
_F_WFOLD = 0.0    # wing/rotor fold weight fraction
_F_MAR   = 0.0    # marinization weight fraction
_F_PRESS = 0.0    # pressurization weight fraction
_F_CW    = 0.06   # crashworthiness weight fraction (of fuselage weight)
_NZ      = 3.5    # design ultimate load factor

_C = -0.0866   # wetted-area regression coefficient
_D = 0.8099    # wetted-area regression exponent coefficient


def compute_fuselage_weight(vehicle_mtow, fuselage_length, tech_factor=1.0,
                             landing_gear_on_fuselage=False, landing_gear_retractable=False,
                             has_cargo_ramp=False):
    """ Calculates fuselage/airframe mass (basic structure, tail/wing fold, marinization,
        pressurization, crashworthiness) using the AFDD84 universal fuselage weight model.

        Source:
            Hydra `afdd/fuselage.py::fuselage_weight`, AFDD84 model. Location/ramp premium factors
            corrected against real NDARC (`weight_model.f90::GetWeightAFDD_Fuselage`) — see module
            docstring.

        Inputs:
            vehicle_mtow                vehicle max takeoff weight    [kg]
            fuselage_length              fuselage length                [m]
            tech_factor                  technology weight-scaling factor    [Unitless]
            landing_gear_on_fuselage     True if the landing gear is fuselage-mounted (activates
                                         the 1.1627 structural premium)    [bool]
            landing_gear_retractable     True if the (fuselage-mounted) gear retracts (activates
                                         the additional 1.1437 premium; no effect if
                                         `landing_gear_on_fuselage` is False, matching NDARC)  [bool]
            has_cargo_ramp               True if the vehicle has a cargo ramp (activates the
                                         1.2749 premium)    [bool]

        Outputs:
            weight:   dict with 'basic', 'tail_folding', 'marinization', 'wing_folding',
                      'pressurization', 'crashworth', 'total', all    [kg]
    """
    gtow = vehicle_mtow * _KG2LB
    l_fus = fuselage_length * _M2F

    f_lgloc = _F_LGLOC_ON_FUSELAGE if landing_gear_on_fuselage else 1.0
    f_lgret = _F_LGRET_RETRACTABLE if (landing_gear_on_fuselage and landing_gear_retractable) else 1.0
    f_ramp = _F_RAMP_CARGO if has_cargo_ramp else 1.0

    s_body = 10**(_C + _D*np.log10(gtow))

    wght_basic = (25.41 * f_lgloc*f_lgret*f_ramp * (gtow*0.001)**0.4879 *
                  (gtow*_NZ*0.001)**0.2075 * (s_body**0.1676) *
                  (l_fus**0.1512))

    wght_tfold = _F_TFOLD * wght_basic
    wght_wfold = _F_WFOLD * (wght_basic + wght_tfold)
    wght_mar = _F_MAR * wght_basic
    wght_press = _F_PRESS * wght_basic

    wght_cw = _F_CW * (wght_basic + wght_tfold + wght_wfold + wght_mar + wght_press)

    wght_basic = wght_basic * tech_factor
    wght_tfold = wght_tfold * tech_factor
    wght_mar = wght_mar * tech_factor
    wght_wfold = wght_wfold * tech_factor
    wght_press = wght_press * tech_factor
    wght_cw = wght_cw * tech_factor

    total = wght_basic + wght_tfold + wght_wfold + wght_mar + wght_press + wght_cw

    return {
        'basic': wght_basic * Units.lb,
        'tail_folding': wght_tfold * Units.lb,
        'marinization': wght_mar * Units.lb,
        'wing_folding': wght_wfold * Units.lb,
        'pressurization': wght_press * Units.lb,
        'crashworth': wght_cw * Units.lb,
        'total': total * Units.lb,
    }
