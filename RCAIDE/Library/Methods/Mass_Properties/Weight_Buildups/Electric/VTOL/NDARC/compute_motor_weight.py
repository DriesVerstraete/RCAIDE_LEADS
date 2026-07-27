# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/NDARC/compute_motor_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

# RCAIDE imports
from RCAIDE.Framework.Core import Units

# ----------------------------------------------------------------------------------------------------------------------
#  Constants
# ----------------------------------------------------------------------------------------------------------------------
# Four independently-sourced motor weight methods, selected by `method`. Three
# (`ndarc`/`electric_motor_data`/`bird2021`) share the identical functional form —
# `coefficient x torque^exponent` — just different calibration constants from different sources;
# kept as separate `method` values rather than folded into one parameterized branch so the source
# of each calibration stays explicit rather than blurred together.
#
# `method='ndarc'` (default): W(Q), the torque-based branch of real NDARC v1.19's
# `GetWeightNASA_Engine_motor` (weight_model.f90, local path:
# `/Users/dverstraete/Sydney Uni Dropbox/Dries Verstraete/Software/NDARC/1_19/NDARC_v1_19_source/`).
# Genuine calibrated regression, confirmed coefficients. NDARC derives its torque input from
# installed power and RPM internally (`q = Power/Omega`, buried in NDARC's own dual metric/
# imperial internal-unit-conversion bookkeeping); ported here to take torque directly instead,
# since RCAIDE already carries `motor.design_torque` and re-deriving NDARC's internal power->torque
# conversion would risk a real unit bug for no benefit — `q = Power/Omega` is dimensionally exact
# either way, this just skips redoing a division RCAIDE's caller has already effectively done.
#
# `method='hydra'` (opt-in, NOT the default): a single term extracted from Hydra's
# `src/Python/Stage_1/afdd/motors.py::weight` (duplicated verbatim in `Stage_1/motors.py`,
# `github.com/VahanaOpenSource/vtol_sizing`) — `wght_motors = nrotors * 2.278*power**0.6563`.
# Flagged explicitly as lower-confidence than the NDARC branch:
#   - That source function is never called anywhere in the Hydra repo (confirmed by full-repo
#     grep) — genuinely dead code, not just under-exercised.
#   - Its own header disclaims it: "the main AFDD empty weight model is still in the fortran."
#   - No `lb2kg`/unit conversion at the end, unlike every sibling `afdd/` file, and no caller
#     anywhere to infer the intended power unit from. Assumed kW here (typical for this class of
#     regression and consistent with the rest of this Hydra branch) — UNVERIFIED, not confirmed
#     against any real NDARC or Hydra source.
#   - Six alternative coefficient sets are commented out immediately next to this one in the
#     source (0.5269*P^0.8983, 2.2*P^0.661, 0.40199028*P^1.15382503, 0.7907*P^0.8246,
#     0.07529*P^1.226) — reads as an unfinished calibration exploration, not a formula the
#     original author settled on and trusted.
# `method='electric_motor_data'`: `SUAVE/Methods/Weights/Buildups/eVTOL/empty_hydra.py`'s own
# inline motor formula (local path: `/Users/dverstraete/Sydney Uni Dropbox/Dries
# Verstraete/SUAVE/`), `motor_coefficient_torque=0.3288`, `motor_exponent_torque=0.7622`, sourced
# per that file's own comment from "electric_motor_data.xls sheet 2" — a real, independently
# calibrated empirical database fit, unrelated to NDARC's Theory Manual regression despite the
# superficially similar coefficients.
#
# `method='bird2021'`: an alternative calibration commented out immediately next to the
# `electric_motor_data` one in the same `empty_hydra.py` source,
# `motor_coefficient_torque=0.4528`, `motor_exponent_torque=0.7224`, citing "Bird2021" (not yet
# traced to a specific paper/report).
#
# The original Hydra function also bundles cable, cooling, battery-approximation, and gearbox
# weight into the same return value (hybrid-turboshaft-generator-architecture specific) — none of
# that is ported here, only the isolated motor term, to avoid double-counting with wiring/battery/
# drive components already handled elsewhere in this method.
#
# See 01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md ("Open
# action items (motor weight)") for the full reconciliation status across all four sources.

_M2F = 1.0 / Units.ft
_KG2LB = 1.0 / Units.lb
_NM2LBFT = 1.0 / 1.3558179483314004   # N*m -> lbf*ft (force-based, not RCAIDE's mass-based Units.lb)

# NDARC W(Q) coefficients, selected by kind_design (torque-to-weight design point):
#   0: "only high Q/W" design point
#   1: "high Q/W" design point (same formula shape as 0, different coefficients)
#   2: "low Q/W" factor — applies a 2.5606 weight multiplier on top of the kind_design=1 formula
_NDARC_COEFF = {
    0: (0.3928, 0.8587, 1.0),
    1: (0.5382, 0.8129, 1.0),
    2: (0.5382, 0.8129, 2.5606),
}

_ELECTRIC_MOTOR_DATA_COEFF = 0.3288   # electric_motor_data.xls sheet 2, via empty_hydra.py
_ELECTRIC_MOTOR_DATA_EXP = 0.7622

_BIRD2021_COEFF = 0.4528              # Bird2021, via empty_hydra.py (commented alternative)
_BIRD2021_EXP = 0.7224

_HYDRA_COEFF = 2.278
_HYDRA_EXP = 0.6563


def _torque_power_law(design_torque, coefficient, exponent):
    """ Shared `coefficient x torque^exponent` evaluator for every torque-based method
        (`ndarc`/`electric_motor_data`/`bird2021`) — same functional form, different calibration
        constants and sources. Torque converted N*m -> lbf*ft internally since all three source
        calibrations are imperial (matches SUAVE's `nasa_motor.py` convention exactly). """
    q = abs(design_torque) * _NM2LBFT
    return coefficient * (q**exponent) * Units.lb


def compute_motor_weight(design_torque=None, motor_power=None, n_rotor=1, method='ndarc',
                          kind_design=0):
    """ Calculates motor mass using one of four independently-sourced methods.

        Source:
            method='ndarc': NDARC v1.19 `weight_model.f90::GetWeightNASA_Engine_motor`, W(Q) branch.
                             Independently confirmed against SUAVE's `nasa_motor.py`.
            method='electric_motor_data': `empty_hydra.py`'s own inline formula, calibrated from
                             "electric_motor_data.xls sheet 2" — not NDARC-derived.
            method='bird2021': alternative calibration from the same `empty_hydra.py` source,
                             citing "Bird2021".
            method='hydra': Hydra `afdd/motors.py::weight`'s `wght_motors` term only — see module
                             docstring for why this branch is lower-confidence than the other three.

        Inputs:
            design_torque    motor design torque — required for method='ndarc'/
                              'electric_motor_data'/'bird2021'                            [N*m]
            motor_power      power per motor — required for method='hydra', UNVERIFIED
                              unit convention, assumed kW                                [kW]
            n_rotor           number of motors (only used by method='hydra' — the torque-based
                              methods are inherently per-motor, sum externally for a group)  [Unitless]
            method            'ndarc' (default), 'electric_motor_data', 'bird2021', or 'hydra' [str]
            kind_design       NDARC torque-to-weight design point, method='ndarc' only:
                              0 (default) or 1 = high Q/W, 2 = low Q/W (heavier)          [int]

        Outputs:
            mass:   motor mass                                                            [kg]
    """
    if method == 'ndarc':
        if design_torque is None:
            raise ValueError("compute_motor_weight: method='ndarc' requires design_torque")
        if kind_design not in _NDARC_COEFF:
            raise ValueError(f"compute_motor_weight: unrecognized kind_design {kind_design!r}")
        k, x, mult = _NDARC_COEFF[kind_design]
        return _torque_power_law(design_torque, k, x) * mult

    elif method == 'electric_motor_data':
        if design_torque is None:
            raise ValueError("compute_motor_weight: method='electric_motor_data' requires design_torque")
        return _torque_power_law(design_torque, _ELECTRIC_MOTOR_DATA_COEFF, _ELECTRIC_MOTOR_DATA_EXP)

    elif method == 'bird2021':
        if design_torque is None:
            raise ValueError("compute_motor_weight: method='bird2021' requires design_torque")
        return _torque_power_law(design_torque, _BIRD2021_COEFF, _BIRD2021_EXP)

    elif method == 'hydra':
        if motor_power is None:
            raise ValueError("compute_motor_weight: method='hydra' requires motor_power")
        wght_motors = n_rotor * _HYDRA_COEFF * (motor_power**_HYDRA_EXP)
        return wght_motors  # no unit conversion in source — see module docstring caveat

    else:
        raise ValueError(
            f"compute_motor_weight: unrecognized method {method!r} "
            "(expected 'ndarc', 'electric_motor_data', 'bird2021', or 'hydra')"
        )
