# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/Hydra/compute_fuselage_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  Constants — Hydra's "GB/ZL" ellipsoid fuselage model
# ----------------------------------------------------------------------------------------------------------------------
# Direct port of Hydra's `src/Python/Stage_1/afdd/fuselage_v2.py::fuselage_v2`
# (`github.com/VahanaOpenSource/vtol_sizing`), confirmed live for the eVTOL (non-helicopter) branch
# of `empty_weight.py:125`. A different physics model from NDARC's AFDD84 regression (sibling
# folder): an ellipsoid wetted-area skin/bulkhead/canopy model plus a keel structure sized against
# three independent load cases (wing-lift bending, wing-torsion, side-landing).
#
# Source is imperial-flavored in its own header ("ALL UNITS IN IMPERIAL") but the actual body is a
# mix — `gtow` arrives in lb from the caller (converted from kg by the source's own rounded
# `kg2lb=1/0.4536`) and is converted to a force in Newtons via `/2.2*9.81` (also rounded), while
# every other length/density/stress quantity (fuselage width/height, wetted area, keel geometry)
# is already SI throughout. Since this port takes `vehicle_mtow` in kg directly (RCAIDE-native),
# the source's kg->lb->N round trip is unnecessary — computed directly as
# `vehicle_mtow * Units.standard_gravity` (exact gravity, 9.80665) instead. Validated 2026-07-28:
# geometry-only outputs (skin/bulkhead/canopy — no weight dependency) match source bit-for-bit;
# `keel` (the only weight-dependent output) carries a ~0.28% relative gap versus source, fully
# explained by the source's own double-rounded unit chain (`kg->lb` via `1/0.4536`, then `lb->N`
# via `/2.2*9.81` — net factor 9.8305 m/s^2 vs the exact 9.80665 used here). A larger gap than the
# ~2e-5 seen on `NDARC`'s wing/rotor/fuselage ports (those only round ft/lb length/mass units, not
# a doubly-rounded gravity chain) but understood and deliberate, not a porting bug.
#
# Fuselage width/height are fixed at 1.55 m each in the source — not exposed as inputs there, so
# not exposed here either (a genuine source limitation: this model cannot represent a
# narrower/wider fuselage cross-section, only length varies).
#
# The source's own `wing`-object dependency (`wing.groups[i].span`/`.lift_frac`/`.nwings`,
# `wing.ngroups`) is abstracted to two explicit scalar inputs here: the maximum span across wing
# groups, and the minimum per-wing lift fraction across wing groups — the only two derived
# quantities the source actually reads off the wing object, so nothing is lost by not carrying
# RCAIDE's full wing collection into this function.

from RCAIDE.Framework.Core import Units

_WIDTH = 1.55     # fuselage width, m — fixed in source, not exposed
_HEIGHT = 1.55    # fuselage height, m — fixed in source, not exposed

_NG = 3.8      # max maneuver load factor
_NL = 3.5      # landing load factor
_SF = 1.5      # safety factor

_UNI_RHO = 1660.0        # uniaxial ply density, kg/m^3
_UNI_STRESS = 450.0e6    # uniaxial ultimate stress, Pa
_BID_RHO = 1660.0        # BID (bi-directional) ply density, kg/m^3
_BID_STRESS = 275.0e6    # BID ultimate compressive stress, Pa (unused directly — kept for parity with source)
_BID_SHEAR = 47.0e6      # BID ultimate shear stress, Pa
_BID_MIN_THK = 0.00046   # minimum gauge, 3-ply
_BID_BEARING = 400.0e6   # bearing allowable, Pa
_CORE_MIN_THK = 0.0064   # core thickness, m
_CORE_RHO = 52.0         # core density, kg/m^3
_PAINT_THK = 0.00015     # m
_PAINT_RHO = 1800.0      # kg/m^3
_CANOPY_THK = 0.003175   # m
_CANOPY_RHO = 1180.0     # kg/m^3
_STEEL_SHEAR = 500.0e6   # bolt shear strength, Pa

_AREAL_WEIGHT = _BID_MIN_THK * _BID_RHO + _CORE_MIN_THK * _CORE_RHO + _PAINT_THK * _PAINT_RHO


def compute_fuselage_weight(vehicle_mtow, wing_span, min_wing_lift_fraction, fuselage_length,
                             tech_factor=1.0, keel_materials=None):
    """ Calculates fuselage mass (skin, bulkheads, canopy, keel) using Hydra's "GB/ZL" ellipsoid
        model: an ellipsoid wetted-area skin/bulkhead/canopy sizing, plus a keel beam sized
        against three independent load cases — wing-lift bending, wing-torsion, and side-landing.

        Source:
            Hydra `afdd/fuselage_v2.py::fuselage_v2`.

        Inputs:
            vehicle_mtow             vehicle max takeoff weight                          [kg]
            wing_span                maximum span across all wing groups (source:
                                       `max(group.span for group in wing.groups)`)        [m]
            min_wing_lift_fraction   minimum per-wing lift fraction across wing groups
                                      (source: `min(group.lift_frac/group.nwings)`)       [Unitless]
            fuselage_length          fuselage length                                     [m]
            tech_factor              technology weight-scaling factor (source applies an
                                       additional flat +20% fastener allowance on top,
                                       always, not gated by tech_factor)                  [Unitless]
            keel_materials           optional RCAIDE `Data()` (or None) overriding the keel
                                       structural material properties below -- 2026-07-30, added
                                       so a real named/sourced `Solid` material class (e.g.
                                       `AS4_3502_Unidirectional_Carbon_Fiber`) can drive the keel
                                       sizing instead of this module's hardcoded constants. Any
                                       field left unset (or `keel_materials=None` entirely) falls
                                       back to the original hardcoded default -- this preserves
                                       bit-for-bit backward compatibility when unset. Recognized
                                       fields, all optional:
                                         uni_stress   0-deg bending/tension allowable       [Pa]
                                         uni_rho      0-deg ply density                     [kg/m^3]
                                         bid_shear    +-45/bidirectional shear allowable    [Pa]
                                         bid_bearing  bidirectional bearing allowable       [Pa]
                                         bid_rho      bidirectional ply density             [kg/m^3]
                                         steel_shear  bolt shear allowable                  [Pa]
                                       Only the keel (structural) terms are overridable -- skin/
                                       bulkhead/canopy sandwich material is a separate, still-
                                       hardcoded assumption (`_AREAL_WEIGHT`, unchanged, not
                                       sourced/revisited this session).

        Outputs:
            weight:   dict with 'skin', 'keel', 'bulkhead', 'canopy', 'total', all        [kg]
    """
    km = keel_materials
    uni_stress  = getattr(km, 'uni_stress',  None) if km is not None else None
    uni_rho     = getattr(km, 'uni_rho',     None) if km is not None else None
    bid_shear   = getattr(km, 'bid_shear',   None) if km is not None else None
    bid_bearing = getattr(km, 'bid_bearing', None) if km is not None else None
    bid_rho     = getattr(km, 'bid_rho',     None) if km is not None else None
    steel_shear = getattr(km, 'steel_shear', None) if km is not None else None

    uni_stress  = _UNI_STRESS   if uni_stress  is None else uni_stress
    uni_rho     = _UNI_RHO      if uni_rho     is None else uni_rho
    bid_shear   = _BID_SHEAR    if bid_shear   is None else bid_shear
    bid_bearing = _BID_BEARING  if bid_bearing is None else bid_bearing
    bid_rho     = _BID_RHO      if bid_rho     is None else bid_rho
    steel_shear = _STEEL_SHEAR  if steel_shear is None else steel_shear

    weight_n = vehicle_mtow * Units.standard_gravity

    width = _WIDTH
    height = _HEIGHT
    length = fuselage_length
    span = wing_span
    min_lf = min_wing_lift_fraction

    # wetted area (ellipsoid approximation) and skin mass
    s_wet = 4 * np.pi * ((((length * width / 4)**1.6 + (length * height / 4)**1.6 +
                            (width * height / 4)**1.6) / 3))**(1 / 1.6)
    m_skin = s_wet * _AREAL_WEIGHT

    bulkhead_mass = 4 * (np.pi * height * width * 0.25) * _AREAL_WEIGHT

    canopy_mass = s_wet / 10 * _CANOPY_THK * _CANOPY_RHO

    # keel mass due to lift
    l_lift = _NG * weight_n * _SF
    m_bend = l_lift * min_lf * length * (2.0 / 3.0)
    beam_width = width
    beam_height = height
    a_bend = m_bend * beam_height / (4 * uni_stress * (beam_height * 0.5)**2)
    mass_keel = a_bend * length * uni_rho

    # keel mass due to wing torsion
    m_torsion = min_lf * l_lift * span * 0.75
    a_torsion = beam_height * beam_width
    t_torsion = 0.5 * m_torsion / (bid_shear * a_torsion) * (_NG * 0.5)
    mass_keel = mass_keel + 2 * (beam_height + beam_width) * t_torsion * bid_rho

    # keel mass due to side landing
    f_landing = _SF * weight_n * _NL * 0.6403
    a_bolt = f_landing / steel_shear
    d_bolt = 2 * np.sqrt(a_bolt / np.pi)
    t_laminate = f_landing / (d_bolt * bid_bearing)
    v_padup = np.pi * (20 * t_laminate)**2 * t_laminate / 3
    mass_keel = mass_keel + 4 * v_padup * bid_rho

    factor = tech_factor * 1.2   # +20% fastener/misc allowance, always applied (matches source)

    mass_skin = m_skin * factor
    mass_keel = mass_keel * factor
    mass_bulkhead = bulkhead_mass * factor
    mass_canopy = canopy_mass * factor

    total = mass_skin + mass_keel + mass_bulkhead + mass_canopy

    return {
        'skin': mass_skin,
        'keel': mass_keel,
        'bulkhead': mass_bulkhead,
        'canopy': mass_canopy,
        'total': total,
    }
