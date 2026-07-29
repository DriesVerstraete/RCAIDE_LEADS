# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/NDARC/compute_operating_empty_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

# RCAIDE
import RCAIDE
from RCAIDE.Framework.Core import Units, Data
import RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.NDARC as NDARC
import RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.Physics_Based as Vahana
import RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.Common as EVTOL_Common
from RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Conventional.Common import compute_payload_weight

# package imports
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  Compute Operating Empty Weight — NDARC method
# ----------------------------------------------------------------------------------------------------------------------
# Orchestrates every validated NDARC component (see
# 01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md for the
# component-by-component porting/validation record) into the same output schema
# `Physics_Based/Vahana` return, so downstream RCAIDE machinery (CG, MOI, AVL mass file, plots)
# works identically regardless of which method is selected.
#
# **New per-vehicle input convention introduced here** — an `.NDARC` namespace on wings/rotors/
# fuselage/vehicle, mirroring `wing.NDARC.tip_propulsor_tags` (introduced with the wing component).
# None of the three current test vehicles set these yet — wiring a real vehicle through
# `method='NDARC'` is a separate follow-up step, not done as part of building this orchestrator.
# Required/optional fields, by object:
#
#   wing.NDARC.tip_propulsor_tags      list[str], REQUIRED per wing/tail that has tilting
#                                       propulsors mounted on it. Absence (or an empty list) means
#                                       this wing has no tip mass and is NOT sized with NDARC's
#                                       Chappell-Peyran method — it falls back to Vahana's
#                                       beam-bending method instead (see "Wing/tail routing" below;
#                                       this is not optional behavior, it's a correctness
#                                       requirement — Chappell-Peyran silently zeroes out primary
#                                       structure for a zero-tip-mass surface, confirmed 2026-07-27
#                                       on tilt_stopped_rotor_v_tail's V-tail).
#   wing.NDARC.tilt_wing               bool, default False (tiltrotor: only the tip nacelle tilts).
#   wing.NDARC.tip_weight_factor       float, default 1.0 (fWtip).
#   wing.NDARC.tip_weight_increment    float, default 0.0 kg (xWtip).
#
#   rotor.NDARC.nu_blade               float, REQUIRED for any wing-mounted or fuselage-mounted
#                                       rotor — blade flap natural frequency ratio (per-rev). No
#                                       RCAIDE Rotor attribute carries this; NDARC's own formulas
#                                       are highly sensitive to it (exponent ~2.5), so no default is
#                                       assumed — raises clearly if missing rather than guessing.
#   rotor.NDARC.aircraft_id            int, default 2 (tiltrotor) — see vehicle.NDARC.aircraft_id
#                                       for the meaning; settable per-rotor in case a vehicle mixes
#                                       rotor types (e.g. tiltrotor + fixed tail rotor).
#
#   vehicle.NDARC.aircraft_id          int, default 2 (tiltrotor). Hydra/NDARC internal
#                                       vehicle-classification flag — 2 activates tilt-factored
#                                       blade/hub weight and the conversion (tilt) flight-control
#                                       boost mechanism; used as the flight-control-system default
#                                       when rotors don't each specify their own.
#   vehicle.NDARC.motor_kind_design    int, default 0 (high torque-to-weight design point).
#   vehicle.NDARC.landing_gear_on_fuselage/retractable/has_cargo_ramp
#                                       bool, all default False. **Every current vehicle will
#                                       realistically need `landing_gear_on_fuselage=True`** once
#                                       wired up (confirmed 2026-07-27: lift_cruise has genuinely
#                                       fuselage-mounted gear; the other two would use fuselage-
#                                       mounted skids) — defaulting to False here anyway rather
#                                       than silently assuming it, since a real vehicle definition
#                                       should state this explicitly, not inherit a guess.
#
# Wing/tail routing: every `vehicle.wings` entry is checked for `wing.NDARC.tip_propulsor_tags`.
# If present and non-empty, sized with NDARC's Chappell-Peyran method (`compute_wing_tip_mass` +
# `compute_wing_weight_group`). If absent or empty, sized with Vahana's beam-bending method
# (`Physics_Based.compute_wing_weight`) instead — this applies to true wings with no tip-mounted
# propulsor and to tail surfaces (RCAIDE routes `Horizontal_Tail`/`Vertical_Tail` into
# `vehicle.wings` too, confirmed 2026-07-27) equally.
#
# Components with no NDARC-specific model (fall back to Vahana/Common, matching the porting
# inventory's own component matrix): wiring, boom, payload, battery/fuel-cell pass-through,
# thermal management, and — per the 2026-07-27 decision — seats/avionics/ECS (NDARC's own versions
# are input-summation structures defaulting to zero, not calibrated regressions; see the inventory
# doc's "Avionics / seats / ECS / systems" row).


def compute_operating_empty_weight(vehicle, settings=None):
    """ Calculates the empty vehicle mass for an eVTOL-type aircraft using NDARC-derived component
        weight methods where a validated NDARC port exists, falling back to `Vahana`/`Physics_Based`
        methods elsewhere. See module docstring for the full per-component method map and the new
        `.NDARC` namespace input convention this orchestrator introduces.

        Inputs:
            vehicle:     RCAIDE Vehicle Data Structure (with `.NDARC` namespace inputs set — see
                         module docstring)
            settings:    Weights analysis settings (same object `Physics_Based` uses;
                         `settings.miscelleneous_weight_factor` applied identically)

        Outputs:
            output:      Data dictionary, same schema as `Physics_Based.compute_operating_empty_weight`
    """
    diff = 100
    iterations = 0
    MTOW = vehicle.mass_properties.max_takeoff
    tolerance = settings.mtow_convergence_tolerance

    v_ndarc = getattr(vehicle, 'NDARC', Data())
    default_aircraft_id = getattr(v_ndarc, 'aircraft_id', 2)
    motor_kind_design = getattr(v_ndarc, 'motor_kind_design', 0)
    lg_on_fuselage = getattr(v_ndarc, 'landing_gear_on_fuselage', False)
    lg_retractable = getattr(v_ndarc, 'landing_gear_retractable', False)
    has_cargo_ramp = getattr(v_ndarc, 'has_cargo_ramp', False)

    # Relative tolerance (fraction of MTOW), not an absolute mass -- matches Physics_Based's
    # 2026-07-27 fix. See settings.mtow_convergence_tolerance docstring in
    # Framework/Analyses/Weights/Weights.py.
    while abs(diff) > tolerance*MTOW:
        miscelleneous_weight_factor = settings.miscelleneous_weight_factor

        weight = Data()
        weight.battery = 0.0
        weight.fuel_cell = 0.0
        weight.payload = 0.0
        weight.servos = 0.0
        weight.hubs = 0.0
        weight.booms = 0.0
        weight.BRS = 0.0
        weight.motors = 0.0
        weight.rotors = 0.0
        weight.fuselage = 0.0
        weight.wiring = 0.0
        weight.wings = Data()
        weight.wings_total = 0.0
        weight.flight_controls = 0.0
        weight.thermal_management_system = Data()

        # -------------------------------------------------------------------------------
        # Payload
        # -------------------------------------------------------------------------------
        payload = compute_payload_weight(vehicle, W_passenger=70. * Units.kg, W_baggage=0 * Units.lbs)
        weight.seats = vehicle.number_of_passengers * 15. * Units.kg
        weight.passengers = payload.passengers
        weight.avionics = 15. * Units.kg
        weight.landing_gear = 0.0  # NDARC's own landing-gear formula, not Vahana's fixed fraction
        weight.ECS = vehicle.number_of_passengers * 7. * Units.kg

        # -------------------------------------------------------------------------------
        # Rotors + motors: blade/hub via NDARC, motor via NDARC (torque-based, default method)
        # -------------------------------------------------------------------------------
        maxLiftPower = 0
        total_number_of_rotors = 0
        eta = 0
        rTip_ref = 0.0
        maxVTip = 0.0
        for network in vehicle.networks:
            for propulsor in network.propulsors:
                rotor = propulsor.rotor
                r_ndarc = getattr(rotor, 'NDARC', Data())
                nu_blade = getattr(r_ndarc, 'nu_blade', None)
                if nu_blade is None:
                    raise ValueError(
                        f"NDARC.compute_operating_empty_weight: rotor '{rotor.tag}' has no "
                        "rotor.NDARC.nu_blade set (blade flap frequency ratio, required — no "
                        "default assumed, see module docstring)."
                    )
                aircraft_id = getattr(r_ndarc, 'aircraft_id', default_aircraft_id)

                n_blade = rotor.number_of_blades
                radius = rotor.tip_radius
                chord = np.mean(rotor.chord_distribution)
                omega = EVTOL_Common.max_design_field(rotor, 'design_angular_velocity')
                tip_speed = omega * radius

                blade_hub, mass_per_unit = NDARC.compute_rotor_weight(
                    aircraft_id=aircraft_id, n_blade=n_blade, n_rotor=1, radius=radius,
                    chord=chord, tip_speed=tip_speed, nu_blade=nu_blade, tech_factor=1.0,
                )
                weight.rotors += blade_hub['blades']
                weight.hubs += blade_hub['hub'] + blade_hub['spinner'] + blade_hub['folding']

                torque = propulsor.motor.design_torque
                motor_mass = NDARC.compute_motor_weight(
                    design_torque=torque, method='empirical_2026', kind_design=motor_kind_design,
                )
                weight.motors += motor_mass
                propulsor.motor.mass_properties.mass = motor_mass
                rotor.mass_properties.mass = blade_hub['total'] if 'total' in blade_hub else (
                    blade_hub['blades'] + blade_hub['hub'] + blade_hub['spinner'] + blade_hub['folding']
                )

                total_number_of_rotors += 1
                eta = propulsor.motor.efficiency
                if radius > rTip_ref:
                    rTip_ref = radius
                    maxVTip = tip_speed

        weight.BRS = NDARC.compute_emergency_system_weight(MTOW, tech_factor=1.0)

        lg = NDARC.compute_landing_gear_weight(MTOW, tech_factor=1.0)
        weight.landing_gear = lg['total']

        # -------------------------------------------------------------------------------
        # Flight control system (rotary-wing + fixed-wing + tilting hardware + hydraulics)
        # -------------------------------------------------------------------------------
        n_wing_surfaces = len(vehicle.wings)
        n_blade_ref = 0
        chord_ref = 0.0
        for network in vehicle.networks:
            for propulsor in network.propulsors:
                n_blade_ref = propulsor.rotor.number_of_blades
                chord_ref = np.mean(propulsor.rotor.chord_distribution)
                break
            break
        fc = NDARC.compute_flight_control_system_weight(
            n_rotor=total_number_of_rotors, n_blade=n_blade_ref, chord=chord_ref,
            tip_speed=maxVTip, vehicle_mtow=MTOW, aircraft_id=default_aircraft_id,
            n_prop=0, n_wing=n_wing_surfaces, rotor_radius=rTip_ref, tech_factor=1.0,
        )
        weight.flight_controls = fc['total']

        # -------------------------------------------------------------------------------
        # Thermal management (pass-through, same as Vahana — no NDARC model)
        # -------------------------------------------------------------------------------
        tms_weight = 0.0
        for network in vehicle.networks:
            for coolant_line in network.coolant_lines:
                for i, battery_module in enumerate(coolant_line.battery_modules):
                    for HAS in battery_module:
                        tms_weight += HAS.mass_properties.mass
                for tag, item in coolant_line.items():
                    if tag == 'heat_exchangers':
                        for heat_exchanger in item:
                            tms_weight += heat_exchanger.mass_properties.mass
                    if tag == 'reservoirs':
                        for reservoir in item:
                            tms_weight += reservoir.mass_properties.mass
        weight.thermal_management_system.total = tms_weight

        # -------------------------------------------------------------------------------
        # Battery / fuel cell (pass-through, same as Vahana)
        # -------------------------------------------------------------------------------
        for network in vehicle.networks:
            for bus in network.busses:
                for modules in bus.battery_modules:
                    weight.battery += modules.mass_properties.mass * Units.kg
                for fuel_cell in bus.fuel_cell_stacks:
                    weight.fuel_cell += fuel_cell.mass_properties.mass * Units.kg

        # -------------------------------------------------------------------------------
        # Wings and tails — NDARC Chappell-Peyran if tip-mass-bearing, NDARC AFDD93 "parametric
        # method" (compute_fixed_wing_weight) otherwise. Both are now genuine NDARC methods, so
        # nothing here falls back to Vahana any more for wings/tails (see 2026-07-28 decision:
        # AFDD93 confirmed source-exact against the real Fortran, replacing the previous Vahana
        # beam-bending fallback for non-tip-mass wings).
        # -------------------------------------------------------------------------------
        maxSpan = 0
        for wing in vehicle.wings:
            maxSpan = max(wing.spans.projected, maxSpan)
            w_ndarc = getattr(wing, 'NDARC', Data())
            tip_tags = getattr(w_ndarc, 'tip_propulsor_tags', [])

            if wing.symbolic:
                wing_weight = 0
            elif tip_tags:
                tip_mass = NDARC.compute_wing_tip_mass(wing, vehicle)
                tag0 = tip_tags[0]
                ref_rotor = None
                for network in vehicle.networks:
                    for propulsor in network.propulsors:
                        if propulsor.tag == tag0:
                            ref_rotor = propulsor.rotor
                if ref_rotor is None:
                    raise ValueError(
                        f"NDARC.compute_operating_empty_weight: wing '{wing.tag}' "
                        f"NDARC.tip_propulsor_tags references '{tag0}', not found on vehicle."
                    )
                tilt_wing = getattr(w_ndarc, 'tilt_wing', False)
                group = NDARC.compute_wing_weight_group(
                    span=wing.spans.projected, area=wing.areas.reference,
                    chord=wing.areas.reference / wing.spans.projected,
                    rotor_angular_velocity=EVTOL_Common.max_design_field(ref_rotor, 'design_angular_velocity'),
                    rotor_radius=ref_rotor.tip_radius, wing_tip_mass=tip_mass,
                    vehicle_mtow=MTOW, tilt_wing=tilt_wing, tech_factor=1.0, n_wings=1,
                )
                wing_weight = group['total']
            else:
                lift_fraction = getattr(w_ndarc, 'lift_fraction', 1.0)
                landing_gear_on_wing = getattr(w_ndarc, 'landing_gear_on_wing', False)
                fixed_wing = NDARC.compute_fixed_wing_weight(
                    vehicle_mtow=MTOW, area=wing.areas.reference, aspect_ratio=wing.aspect_ratio,
                    taper=wing.taper, thickness_to_chord=wing.thickness_to_chord,
                    lift_fraction=lift_fraction, sweep=wing.sweeps.quarter_chord,
                    landing_gear_on_wing=landing_gear_on_wing, tech_factor=1.0,
                )
                wing_weight = fixed_wing['total']

            wing_tag = wing.tag
            weight.wings[wing_tag] = wing_weight
            wing.mass_properties.mass = wing_weight
            weight.wings_total += wing_weight

            if isinstance(wing, RCAIDE.Library.Components.Wings.Main_Wing):
                wiring_weight = EVTOL_Common.compute_wiring_weight(
                    wing, vehicle, maxLiftPower / (eta * total_number_of_rotors) if eta and total_number_of_rotors else 0.0
                ) * Units.kg
            else:
                wiring_weight = 0
            weight.wiring += wiring_weight

        # -------------------------------------------------------------------------------
        # Fuselage — NDARC AFDD84, with real landing-gear-location/ramp premiums
        # -------------------------------------------------------------------------------
        for fuse in vehicle.fuselages:
            fus = NDARC.compute_fuselage_weight(
                vehicle_mtow=MTOW, fuselage_length=fuse.lengths.total, tech_factor=1.0,
                landing_gear_on_fuselage=lg_on_fuselage, landing_gear_retractable=lg_retractable,
                has_cargo_ramp=has_cargo_ramp,
            )
            fuse.mass_properties.center_of_gravity[0][0] = .45 * fuse.lengths.total
            fuse.mass_properties.mass = fus['total']
            weight.fuselage += fus['total']

        # -------------------------------------------------------------------------------
        # Boom (no NDARC model — Vahana/RCAIDE only)
        # -------------------------------------------------------------------------------
        for boom in vehicle.booms:
            boom_weight = Vahana.compute_boom_weight(boom) * Units.kg
            weight.booms += boom_weight
            boom.mass_properties.mass = boom_weight

        # -------------------------------------------------------------------------------
        # Pack up outputs — same schema as Physics_Based
        # -------------------------------------------------------------------------------
        output = Data()
        output.empty = Data()

        output.empty.structural = Data()
        output.empty.structural.wings = miscelleneous_weight_factor * weight.wings_total
        output.empty.structural.fuselage = miscelleneous_weight_factor * weight.fuselage
        output.empty.structural.empennage = 0.0
        output.empty.structural.landing_gear = miscelleneous_weight_factor * weight.landing_gear
        output.empty.structural.nacelle = 0.0
        output.empty.structural.booms = miscelleneous_weight_factor * weight.booms
        output.empty.structural.paint = 0.0
        output.empty.structural.total = (output.empty.structural.wings + output.empty.structural.fuselage
                                          + output.empty.structural.empennage + output.empty.structural.landing_gear
                                          + output.empty.structural.nacelle + output.empty.structural.booms
                                          + output.empty.structural.paint)

        output.empty.propulsion = Data()
        output.empty.propulsion.engines = miscelleneous_weight_factor * weight.rotors
        output.empty.propulsion.thrust_reversers = 0.0
        output.empty.propulsion.miscellaneous = miscelleneous_weight_factor * (weight.BRS + weight.fuel_cell)
        output.empty.propulsion.fuel_system = 0.0
        output.empty.propulsion.fuel_tanks = 0.0
        output.empty.propulsion.electrical_cabling = miscelleneous_weight_factor * weight.wiring
        output.empty.propulsion.thermal_management = miscelleneous_weight_factor * weight.thermal_management_system.total
        output.empty.propulsion.battery = miscelleneous_weight_factor * weight.battery
        output.empty.propulsion.motors = miscelleneous_weight_factor * (weight.motors + weight.servos + weight.hubs)
        output.empty.propulsion.total = (output.empty.propulsion.engines + output.empty.propulsion.thrust_reversers
                                          + output.empty.propulsion.miscellaneous + output.empty.propulsion.fuel_system
                                          + output.empty.propulsion.electrical_cabling + output.empty.propulsion.thermal_management
                                          + output.empty.propulsion.battery + output.empty.propulsion.motors)

        output.empty.systems = Data()
        output.empty.systems.control_systems = miscelleneous_weight_factor * weight.flight_controls
        output.empty.systems.apu = 0.0
        output.empty.systems.electrical = 0.0
        output.empty.systems.avionics = miscelleneous_weight_factor * weight.avionics
        output.empty.systems.hydraulics = 0.0
        output.empty.systems.furnishings = miscelleneous_weight_factor * weight.seats
        output.empty.systems.air_conditioner = miscelleneous_weight_factor * weight.ECS
        output.empty.systems.instruments = 0.0
        output.empty.systems.total = (output.empty.systems.control_systems + output.empty.systems.apu
                                       + output.empty.systems.electrical + output.empty.systems.avionics
                                       + output.empty.systems.hydraulics + output.empty.systems.furnishings
                                       + output.empty.systems.air_conditioner + output.empty.systems.instruments)

        output.payload = Data()
        output.payload.passengers = weight.passengers
        output.payload.baggage = 0.0
        output.payload.cargo = weight.payload
        output.payload.total = output.payload.passengers + output.payload.baggage + output.payload.cargo

        output.operational_items = Data()
        output.operational_items.misc = 0.0
        output.operational_items.flight_crew = 0.0
        output.operational_items.flight_attendants = 0.0
        output.operational_items.passenger_service = 0.0
        output.operational_items.total = (output.operational_items.misc + output.operational_items.flight_crew
                                           + output.operational_items.flight_attendants + output.operational_items.passenger_service)

        output.empty.total = output.empty.systems.total + output.empty.propulsion.total + output.empty.structural.total + output.operational_items.total
        output.zero_fuel_weight = output.empty.total + output.payload.total
        output.max_takeoff = output.empty.total + output.payload.total

        diff = MTOW - output.max_takeoff
        MTOW -= diff
        iterations += 1

        if iterations == 100:
            print('Weight convergence failed!')
            return output

    vehicle.mass_properties.max_takeoff = output.max_takeoff
    vehicle.mass_properties.takeoff = output.max_takeoff

    return output
