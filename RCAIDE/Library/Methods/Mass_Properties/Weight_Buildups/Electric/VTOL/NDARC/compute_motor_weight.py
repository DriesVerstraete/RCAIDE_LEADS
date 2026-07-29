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
# Five independently-sourced motor weight methods, selected by `method`. Four
# (`empirical_2026`/`ndarc`/`spl`/`bird2021`) share the identical functional form —
# `coefficient x torque^exponent` — just different calibration constants from different sources;
# kept as separate `method` values rather than folded into one parameterized branch so the source
# of each calibration stays explicit rather than blurred together. `empirical_2026` is the default
# as of 2026-07-29 (see its own entry below) — `ndarc` remains available and is not deprecated.
#
# `method='ndarc'`: W(Q), the torque-based branch of real NDARC v1.19's
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
# `method='spl'`: `motor_coefficient_torque=0.3288`, `motor_exponent_torque=0.7622` — a real,
# independently-calibrated empirical database fit, unrelated to NDARC's Theory Manual regression
# despite the superficially similar coefficients. Sourced from the lab's own motor database,
# `electric_motor_data.xlsx` sheet 2 (local copies:
# `/Users/dverstraete/Sydney Uni Dropbox/Dries Verstraete/SUAVE/electric_motors/electric_motor_data.xlsx`,
# also `PythonOptimisation/eVTOLAircraftDesign/validation_data/electric_motor_data.xlsx`). These
# exact coefficients also appear, already fit, inline in `SUAVE/Methods/Weights/Buildups/eVTOL/
# empty_hydra.py` (local path: `/Users/dverstraete/Sydney Uni Dropbox/Dries Verstraete/SUAVE/`) —
# that file is where this port's coefficients were taken from directly, since the fit itself isn't
# re-derived from the raw spreadsheet here.
#
# `method='bird2021'`: an alternative calibration commented out immediately next to the `spl` one
# in the same `empty_hydra.py` source, `motor_coefficient_torque=0.4528`,
# `motor_exponent_torque=0.7224`, citing "Bird2021" (not yet traced to a specific paper/report).
#
# The original Hydra function also bundles cable, cooling, battery-approximation, and gearbox
# weight into the same return value (hybrid-turboshaft-generator-architecture specific) — none of
# that is ported here, only the isolated motor term, to avoid double-counting with wiring/battery/
# drive components already handled elsewhere in this method.
#
# See 01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md ("Open
# action items (motor weight)") for the full reconciliation status across all four sources.
#
# `method='empirical_2026'` (DEFAULT as of 2026-07-29, supersedes `'ndarc'` as default): a fresh
# log-log power-law regression fit directly against the lab's own `electric_motor_data.xlsx`
# (`Sydney Uni Dropbox/Dries Verstraete/SUAVE/electric_motors/electric_motor_data.xlsx`, `Sheet1`),
# the same spreadsheet `'spl'`/`'bird2021'` were originally derived from — but fit fresh this
# session (2026-07-29) against 30 real motors' CONTINUOUS torque and weight (not max/peak torque,
# and not the smaller subset `'spl'`/`'bird2021'` used), plus two new real eVTOL propulsion motors
# not in the spreadsheet (BETA Technologies H500B/V600A, `beta.team/motor`). One outlier excluded
# (Evo motor AF340 — real motor, genuinely far lighter than its torque class vs. every other motor
# in the set; flagged not deleted, worth an independent datasheet check). Fit: `mass = 0.2324 *
# Q^0.8235` (Q already in N*m, mass directly in kg — **no Nm->lbf-ft conversion**, unlike
# `ndarc`/`spl`/`bird2021` below, since this fit was performed directly in SI units against the
# spreadsheet's native columns). n=29, R^2=0.9774, mean|error|=19.0% on the fit's own training set
# — markedly better than `'ndarc'` (mean|err|=33.0%) or `'spl'` (mean|err|=53.1%) checked against
# the same real dataset. A parallel power-based fit (`mass = 0.3540 * P^0.9559`, continuous power)
# was also tried and rejected: needed 5 of 29 points excluded (vs. this fit's 1 of 30) to reach
# comparable R^2, and per-motor errors between the two fits disagree by a mean 38 percentage
# points even where aggregate errors look similar — torque is the more robust single predictor.
# Averaging the two fits was also tested and found to be WORSE than this fit alone (contaminated
# by the power fit's occasional large misses) — do not blend torque- and power-based predictions.
# Full derivation, outlier list, and the rejected power-based/blended alternatives:
# 01-mission-profiles/00-decisions/2026-07-29-rotor-motor-weight-formula-comparison.md.

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

_SPL_COEFF = 0.3288   # SPL's electric_motor_data.xlsx sheet 2, via empty_hydra.py
_SPL_EXP = 0.7622

_BIRD2021_COEFF = 0.4528              # Bird2021, via empty_hydra.py (commented alternative)
_BIRD2021_EXP = 0.7224

_HYDRA_COEFF = 2.278
_HYDRA_EXP = 0.6563

_EMPIRICAL_2026_COEFF = 0.2324   # fit directly in SI (Nm -> kg), NOT imperial like the three above
_EMPIRICAL_2026_EXP = 0.8235


def _torque_power_law(design_torque, coefficient, exponent):
    """ Shared `coefficient x torque^exponent` evaluator for every torque-based method
        (`ndarc`/`spl`/`bird2021`) — same functional form, different calibration
        constants and sources. Torque converted N*m -> lbf*ft internally since all three source
        calibrations are imperial (matches SUAVE's `nasa_motor.py` convention exactly). """
    q = abs(design_torque) * _NM2LBFT
    return coefficient * (q**exponent) * Units.lb


def compute_motor_weight(design_torque=None, motor_power=None, n_rotor=1, method='empirical_2026',
                          kind_design=0):
    """ Calculates motor mass using one of five independently-sourced methods.

        Source:
            method='empirical_2026' (DEFAULT): fresh log-log power-law regression fit directly
                             against the lab's own `electric_motor_data.xlsx` plus two new real
                             eVTOL motor points (BETA H500B/V600A) — see module docstring for
                             full derivation, dataset, and why it was chosen over `'ndarc'`.
            method='ndarc': NDARC v1.19 `weight_model.f90::GetWeightNASA_Engine_motor`, W(Q) branch.
                             Independently confirmed against SUAVE's `nasa_motor.py`.
            method='spl': SPL's own `electric_motor_data.xlsx` sheet 2 calibration (coefficients
                             taken from `empty_hydra.py`'s inline copy of the fit) — not NDARC-derived.
            method='bird2021': alternative calibration from the same `empty_hydra.py` source,
                             citing "Bird2021".
            method='hydra': Hydra `afdd/motors.py::weight`'s `wght_motors` term only — see module
                             docstring for why this branch is lower-confidence than the other three.

        Inputs:
            design_torque    motor design torque — required for method='empirical_2026'/'ndarc'/
                              'spl'/'bird2021'                                            [N*m]
            motor_power      power per motor — required for method='hydra', UNVERIFIED
                              unit convention, assumed kW                                [kW]
            n_rotor           number of motors (only used by method='hydra' — the torque-based
                              methods are inherently per-motor, sum externally for a group)  [Unitless]
            method            'empirical_2026' (default), 'ndarc', 'spl', 'bird2021', or 'hydra' [str]
            kind_design       NDARC torque-to-weight design point, method='ndarc' only:
                              0 (default) or 1 = high Q/W, 2 = low Q/W (heavier)          [int]

        Outputs:
            mass:   motor mass                                                            [kg]
    """
    if method == 'empirical_2026':
        if design_torque is None:
            raise ValueError("compute_motor_weight: method='empirical_2026' requires design_torque")
        q = abs(design_torque)
        return _EMPIRICAL_2026_COEFF * (q ** _EMPIRICAL_2026_EXP)

    elif method == 'ndarc':
        if design_torque is None:
            raise ValueError("compute_motor_weight: method='ndarc' requires design_torque")
        if kind_design not in _NDARC_COEFF:
            raise ValueError(f"compute_motor_weight: unrecognized kind_design {kind_design!r}")
        k, x, mult = _NDARC_COEFF[kind_design]
        return _torque_power_law(design_torque, k, x) * mult

    elif method == 'spl':
        if design_torque is None:
            raise ValueError("compute_motor_weight: method='spl' requires design_torque")
        return _torque_power_law(design_torque, _SPL_COEFF, _SPL_EXP)

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
            "(expected 'empirical_2026', 'ndarc', 'spl', 'bird2021', or 'hydra')"
        )
