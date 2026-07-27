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
# is a faithful port of the main formula. Two real discrepancies found and ported faithfully
# (matching Hydra, not "corrected" to match NDARC) rather than silently:
#   1. `_F_LGLOC`/`_F_LGRET`/`_F_RAMP` are hardcoded to 1.0 (off) here, same as Hydra — real NDARC
#      switches these based on actual landing-gear placement/retraction and cargo-ramp presence
#      (`place_LG`, `kind_LG`, `KIND_ramp`). Correct for a vehicle with wing-mounted gear and no
#      ramp (matches NDARC's own default in that case) — but there's no way to activate them here
#      if a future vehicle needs fuselage-mounted or retractable gear.
#   2. Fold-weight base quantity does NOT match NDARC. This port (and Hydra) computes tail-fold
#      weight as a fraction of *fuselage* weight and wing-fold weight as a fraction of
#      *(fuselage+tail-fold)* weight. Real NDARC computes tail-fold as a fraction of TAIL weight
#      and wing-fold as a fraction of WING weight — a different physical basis. Currently dormant
#      (`_F_TFOLD`/`_F_WFOLD` are both 0.0), so numerically inert — but wrong if fold weight is
#      ever activated using these fractions. Needs a real fix (passing tail/wing weight in) before
#      fold weight can be trusted, not just a coefficient tweak.

_M2F = 1.0 / Units.ft
_KG2LB = 1.0 / Units.lb

_F_LGLOC = 1.0    # landing-gear-on-fuselage location factor
_F_LGRET = 1.0    # retractable landing gear factor
_F_RAMP  = 1.0    # cargo ramp factor (1.0 = no ramp)
_F_TFOLD = 0.0    # tail fold weight fraction
_F_WFOLD = 0.0    # wing/rotor fold weight fraction
_F_MAR   = 0.0    # marinization weight fraction
_F_PRESS = 0.0    # pressurization weight fraction
_F_CW    = 0.06   # crashworthiness weight fraction (of fuselage weight)
_NZ      = 3.5    # design ultimate load factor

_C = -0.0866   # wetted-area regression coefficient
_D = 0.8099    # wetted-area regression exponent coefficient


def compute_fuselage_weight(vehicle_mtow, fuselage_length, tech_factor=1.0):
    """ Calculates fuselage/airframe mass (basic structure, tail/wing fold, marinization,
        pressurization, crashworthiness) using the AFDD84 universal fuselage weight model.

        Source:
            Hydra `afdd/fuselage.py::fuselage_weight`, AFDD84 model.

        Inputs:
            vehicle_mtow       vehicle max takeoff weight    [kg]
            fuselage_length    fuselage length                [m]
            tech_factor        technology weight-scaling factor    [Unitless]

        Outputs:
            weight:   dict with 'basic', 'tail_folding', 'marinization', 'wing_folding',
                      'pressurization', 'crashworth', 'total', all    [kg]
    """
    gtow = vehicle_mtow * _KG2LB
    l_fus = fuselage_length * _M2F

    s_body = 10**(_C + _D*np.log10(gtow))

    wght_basic = (25.41 * _F_LGLOC*_F_LGRET*_F_RAMP * (gtow*0.001)**0.4879 *
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
