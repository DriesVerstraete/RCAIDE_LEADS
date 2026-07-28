# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/Hydra/compute_motor_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  Compute motor weight
# ----------------------------------------------------------------------------------------------------------------------
# Direct port of Hydra's `src/Python/Stage_1/engines/DC_motor.py::DC_motor`
# (`github.com/VahanaOpenSource/vtol_sizing`) — NOT `afdd/motors.py::weight`
# (that function is confirmed dead code, never called anywhere in the Hydra repo; it is ported
# separately, opt-in only, as `method='hydra'` in `NDARC/compute_motor_weight.py`).
#
# `DC_motor` is Hydra's actual, live Stage_1 motor-weight path: traced via the real orchestrator
# call chain `empty_weight.py` -> `transmission.weight_rollup()` (`transmission_class.py:31`,
# confirmed the live weight-buildup entry point, not the module-level `transmissions.py` functions
# of the same name) -> `transmission_group.weight()` (`type=='electric'` branch) ->
# `motor_sizing()` (`transmission_class.py:176`) -> `DC_motor(P_req) * tech_factor`, then
# multiplied by `nmotors` by the caller. Confirmed 2026-07-28 by grepping every caller of
# `DC_motor`/`electric_motor`/`weight_rollup` in the Hydra repo — this is the only path that is
# both non-dead and reaches the real empty-weight rollup (`electric_motor.getWeight()`, the other
# live caller of `DC_motor`, is used only in `Postprocessing/vehicle_performance.py` — post-hoc
# performance evaluation, not the Stage_1 weight buildup).

def compute_motor_weight(motor_power, tech_factor=1.0, n_rotor=1):
    """ Calculates motor mass (including electronic speed controller) from rated power per motor,
        using Hydra's piecewise `DC_motor` fit.

        Source:
            Hydra `engines/DC_motor.py::DC_motor`.
            P <= 10 kW: m_motor = 0.2513 x P (small-scale fit), plus a separate ESC term
                        m_esc = (P/12)x32/7.25 (current-based, from an 8 lb QBiT design point).
                        Motor and ESC are summed; source does not model them as separable outputs.
            P > 10 kW:  m_total = 0.8373 x P^0.726 (power-law fit; ESC weight already folded in).

        Inputs:
            motor_power    rated power per motor    [kW]
            tech_factor    technology weight-scaling factor, applied to the per-motor mass
                            (matches source: `motor_sizing()` applies `tech_factor` after
                            `DC_motor()`, before multiplying by `nmotors`)    [Unitless]
            n_rotor        number of motors (source: `transmission_class.py` multiplies the
                            per-motor mass by `self.nmotors` after `tech_factor`)    [Unitless]

        Outputs:
            mass:          total motor (+ ESC) mass for all `n_rotor` motors    [kg]
    """
    if motor_power <= 10.0:
        m_motor = 0.2513 * motor_power
        current = motor_power / 12.0     # kilo-amps, at an assumed 12 V setting
        m_esc = current * 32.0 / 7.25
        m_total = m_motor + m_esc
    else:
        m_total = 0.8373 * motor_power**0.726

    return m_total * tech_factor * n_rotor
