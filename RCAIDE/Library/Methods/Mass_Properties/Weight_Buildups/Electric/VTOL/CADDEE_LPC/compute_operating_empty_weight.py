# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/CADDEE_LPC/compute_operating_empty_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

import RCAIDE
from RCAIDE.Framework.Core import Units, Data
import RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.Physics_Based as EVTOL
import RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.CADDEE_LPC as CADDEE_LPC
import RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.Common as EVTOL_Common
from RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Conventional.Common import compute_payload_weight

import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  Compute Operating Empty Weight — CADDEE_LPC method
# ----------------------------------------------------------------------------------------------------------------------
# Direct copy of `Physics_Based/compute_operating_empty_weight.py` — same pattern already used for
# `NDARC`/`Hydra` (both are copies with component calls swapped in-place, not calls-out to a
# separate decoupled pass). First attempt at this file (2026-07-28) ran a single, one-off
# `Vahana` pass to get rotor/motor/BRS/landing-gear/etc, then overlaid `CADDEE_LPC`'s structural
# numbers on top — WRONG, caught by the user before being finalized: landing gear (2% of MTOW) and
# rotor structural mass (scales with `maxLift ∝ MTOW`) would then be computed at Vahana's own
# (lower) converged MTOW rather than the true, heavier, self-consistent MTOW under CADDEE_LPC's
# structural assumptions — confirmed by direct measurement on `lift_cruise`: Vahana-only converges
# to 2613.8kg vs. CADDEE_LPC's true 2773.7kg (~6% higher), so a decoupled pass would have quietly
# under-sized landing gear and rotor structure. Fixed by keeping everything in the same converging
# loop instead, exactly like `NDARC`/`Hydra` do.
#
# Only wing/fuselage/boom (`EVTOL.compute_*` -> `CADDEE_LPC.compute_*`) are swapped; empennage is
# a new addition Vahana doesn't have its own line for (its wing loop just runs h-tail/v-tail
# through the same wing function as the main wing — CADDEE_LPC computes them from a dedicated
# regression, added as a new bucket, not replacing anything). Rotor/motor/BRS/landing-gear/wiring/
# thermal-management/battery/systems are untouched from `Physics_Based` — CADDEE_alpha has no
# model for any of those (confirmed by direct repo read, 2026-07-28; see
# `20-rcaide-weight-method-porting-inventory.md`).
#
# **This method is valid ONLY for lift+cruise-class vehicles** — a fixed main wing (no tip-mounted
# tilting rotor), boom-mounted lift rotors, a single fuselage, and h-tail+v-tail empennage,
# resembling the NASA LPC reference vehicle the regression was fit to. Nothing here enforces that
# at runtime — it is the caller's responsibility, same as choosing a sensible `aircraft_id` for
# `NDARC`.
#
# **New per-vehicle input convention introduced here** — a `.CADDEE_LPC` namespace on the vehicle,
# for the 2 design variables the nasa_lpc regression needs that don't map onto any single RCAIDE
# component (`wing_area`/`wing_AR`/`fuselage_length` are read directly off the real wing/fuselage
# components instead, no namespace needed for those):
#
#   vehicle.CADDEE_LPC.battery_mass       float, REQUIRED — kg. Not read from the vehicle's own
#                                          battery modules automatically, since the regression was
#                                          fit against a single scalar design variable, not a sum
#                                          over an arbitrary battery configuration. No default.
#   vehicle.CADDEE_LPC.cruise_speed       float, REQUIRED — m/s, the vehicle's design cruise speed.
#                                          No default assumed (same reasoning as battery_mass).
#   vehicle.CADDEE_LPC.correction_factor  float, default 1.3 — applied to ALL FOUR components
#                                          (wing/fuselage/boom/empennage) uniformly. CADDEE_alpha's
#                                          own reference case only ever validates this factor for
#                                          fuselage/boom/empennage (never wing, which they compute
#                                          via real `aframe` FEA instead) — see the inventory doc's
#                                          numeric cross-check before trusting the wing number at
#                                          this default.


def compute_operating_empty_weight(vehicle, settings=None):
    """ Calculates the empty vehicle mass for a lift+cruise-class eVTOL aircraft, using
        LSDOlab/UCSD's `CADDEE_alpha` nasa_lpc regression for wing/fuselage/boom/empennage mass,
        and `Vahana`/`Physics_Based`'s own methods for everything else (seats, avionics,
        servomotors, ballistic recovery system, rotor and hub assembly, landing gear) —
        CADDEE_alpha has no model for any of those. See module docstring for the full method map,
        the `.CADDEE_LPC` namespace convention, and the lift+cruise-only scope.

        Inputs:
            vehicle:     RCAIDE Vehicle Data Structure (with `.CADDEE_LPC` namespace inputs set —
                         see module docstring). Must be a lift+cruise-class layout.
            settings:    Weights analysis settings (same object every other method uses;
                         `settings.miscelleneous_weight_factor` applied identically)

        Outputs:
            output:      Data dictionary, same schema as `Physics_Based.compute_operating_empty_weight`
    """
    v_caddee = getattr(vehicle, 'CADDEE_LPC', Data())
    battery_mass_design_var = getattr(v_caddee, 'battery_mass', None)
    cruise_speed = getattr(v_caddee, 'cruise_speed', None)
    correction_factor = getattr(v_caddee, 'correction_factor', 1.3)

    if battery_mass_design_var is None:
        raise ValueError(
            "CADDEE_LPC.compute_operating_empty_weight: vehicle.CADDEE_LPC.battery_mass not set "
            "(required — no default assumed, see module docstring)."
        )
    if cruise_speed is None:
        raise ValueError(
            "CADDEE_LPC.compute_operating_empty_weight: vehicle.CADDEE_LPC.cruise_speed not set "
            "(required — no default assumed, see module docstring)."
        )

    main_wing = None
    for wing in vehicle.wings:
        if isinstance(wing, RCAIDE.Library.Components.Wings.Main_Wing):
            main_wing = wing
    if main_wing is None:
        raise ValueError(
            "CADDEE_LPC.compute_operating_empty_weight: no Main_Wing found on vehicle — this "
            "method requires a lift+cruise-class layout (see module docstring)."
        )

    fuselage_length_for_regression = None
    for fuse in vehicle.fuselages:
        fuselage_length_for_regression = fuse.lengths.total
        break
    if fuselage_length_for_regression is None:
        raise ValueError(
            "CADDEE_LPC.compute_operating_empty_weight: no fuselage found on vehicle."
        )

    wing_area_for_regression = main_wing.areas.reference
    wing_AR_for_regression = main_wing.aspect_ratio

    diff = 100
    iterations = 0
    MTOW = vehicle.mass_properties.max_takeoff
    tolerance = settings.mtow_convergence_tolerance

    while abs(diff) > tolerance*MTOW:

        miscelleneous_weight_factor = settings.miscelleneous_weight_factor
        safety_factor = 1.5
        disk_area_factor = 1.15
        max_thrust_to_weight_ratio = 1.1
        max_g_load = 3.8

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
        weight.empennage_total = 0.0
        weight.thermal_management_system = Data()

        atmosphere = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
        atmo_data = atmosphere.compute_values(0, 0)
        rho_ref = atmo_data.density[0, 0]
        maxLift = MTOW * max_thrust_to_weight_ratio * 9.81
        AvgBladeCD = 0.012

        # -------------------------------------------------------------------------------
        # Payload
        # -------------------------------------------------------------------------------
        payload = compute_payload_weight(vehicle, W_passenger=70. * Units.kg, W_baggage=0 * Units.lbs)

        weight.seats = vehicle.number_of_passengers * 15. * Units.kg
        weight.passengers = payload.passengers
        weight.avionics = 15. * Units.kg
        weight.landing_gear = MTOW * 0.02 * Units.kg
        weight.ECS = vehicle.number_of_passengers * 7. * Units.kg

        # -------------------------------------------------------------------------------
        # Network Weight — rotor/motor/BRS, untouched from Physics_Based
        # -------------------------------------------------------------------------------
        maxLiftPower = 0
        total_number_of_rotors = 0
        maxVTip = 0
        eta = 0
        for network in vehicle.networks:
            for system in network.systems:
                if type(system) == RCAIDE.Library.Components.Powertrain.Systems.Avionics:
                    weight.avionics += system.mass_properties.mass * Units.kg
                if type(system) == RCAIDE.Library.Components.Powertrain.Systems.Environmental_Controls:
                    system.mass_properties.mass = weight.ECS
                if type(system) == RCAIDE.Library.Components.Powertrain.Systems.Electrical:
                    system.mass_properties.mass = weight.wiring
                if type(system) == RCAIDE.Library.Components.Powertrain.Systems.Furnishings:
                    system.mass_properties.mass = weight.seats

            for bus in network.busses:
                for modules in bus.battery_modules:
                    weight.battery += modules.mass_properties.mass * Units.kg
                for fuel_cell in bus.fuel_cell_stacks:
                    weight.fuel_cell += fuel_cell.mass_properties.mass * Units.kg

                lift_rotor_hub_weight = 4. * Units.kg
                prop_hub_weight = 4. * Units.kg
                lift_rotor_BRS_weight = 16. * Units.kg

                number_of_propellers = 0.0
                number_of_lift_rotors = 0.0
                total_number_of_rotors = 0.0
                lift_rotor_servo_weight = 0.0

            for propulsor in network.propulsors:
                rotor = propulsor.rotor
                if type(rotor) == RCAIDE.Library.Components.Powertrain.Converters.Propeller:
                    number_of_propellers += 1
                    rTip_ref = rotor.tip_radius
                    bladeSol_ref = rotor.blade_solidity
                    prop_servo_weight = 5.2 * Units.kg
                    propeller_mass = EVTOL_Common.compute_rotor_weight(rotor, maxLift/5.) * Units.kg
                    weight.rotors += propeller_mass
                    rotor.mass_properties.mass = propeller_mass + prop_hub_weight + prop_servo_weight
                    maxVTip = rotor.cruise.design_angular_velocity * rotor.tip_radius
                    weight.servos += prop_servo_weight
                    weight.hubs += prop_hub_weight

                if (type(rotor) == RCAIDE.Library.Components.Powertrain.Converters.Lift_Rotor or
                        type(rotor) == RCAIDE.Library.Components.Powertrain.Converters.Prop_Rotor) or \
                        type(rotor) == RCAIDE.Library.Components.Powertrain.Converters.Rotor:
                    number_of_lift_rotors += 1
                    rTip_ref = rotor.tip_radius
                    bladeSol_ref = rotor.blade_solidity
                    maxVTip = EVTOL_Common.max_design_field(rotor, 'design_angular_velocity') * rotor.tip_radius
                    lift_rotor_servo_weight = 0.65 * Units.kg
                    if rotor.oei.design_thrust is None:
                        design_thrust = EVTOL_Common.max_design_field(rotor, 'design_thrust')
                    else:
                        design_thrust = rotor.oei.design_thrust
                    lift_rotor_mass = EVTOL_Common.compute_rotor_weight(rotor, design_thrust)
                    weight.rotors += lift_rotor_mass
                    rotor.mass_properties.mass = lift_rotor_mass + lift_rotor_hub_weight + lift_rotor_servo_weight
                    weight.servos += lift_rotor_servo_weight
                    weight.hubs += lift_rotor_hub_weight

                eta = propulsor.motor.efficiency
                weight.motors += propulsor.motor.mass_properties.mass

            total_number_of_rotors = int(number_of_lift_rotors + number_of_propellers)
            if total_number_of_rotors > 1:
                prop_BRS_weight = 16. * Units.kg
            else:
                prop_BRS_weight = 0. * Units.kg

            weight.BRS += (prop_BRS_weight + lift_rotor_BRS_weight)
            maxLiftPower = 1.15*maxLift*(disk_area_factor*np.sqrt(maxLift/(2*rho_ref*np.pi*rTip_ref**2)) +
                                          bladeSol_ref*AvgBladeCD/8*maxVTip**3/(maxLift/(rho_ref*np.pi*rTip_ref**2)))

            if number_of_lift_rotors == 1:
                maxLiftOmega = maxVTip/rTip_ref
                maxLiftTorque = maxLiftPower / maxLiftOmega
                for bus in network.busses:
                    tailrotor = next(iter(bus.lift_rotors))
                    weight.tail_rotor = EVTOL_Common.compute_rotor_weight(tailrotor, 1.5*maxLiftTorque/(1.25*rTip_ref))*0.2 * Units.kg
                    weight.rotors += weight.tail_rotor

            tms_weight = 0.0
            for coolant_line in network.coolant_lines:
                weight.thermal_management_system.battery_module = Data()
                for i, battery_module in enumerate(coolant_line.battery_modules):
                    module_key = f'module_{i+1}'
                    weight.thermal_management_system.battery_module[module_key] = 0.0
                    for HAS in battery_module:
                        weight.thermal_management_system.battery_module[module_key] = HAS.mass_properties.mass
                        tms_weight += HAS.mass_properties.mass

                for tag, item in coolant_line.items():
                    if tag == 'heat_exchangers':
                        for heat_exchanger in item:
                            weight.thermal_management_system[heat_exchanger.tag] = heat_exchanger.mass_properties.mass
                            tms_weight += heat_exchanger.mass_properties.mass
                    if tag == 'reservoirs':
                        for reservoir in item:
                            weight.thermal_management_system[reservoir.tag] = reservoir.mass_properties.mass
                            tms_weight += reservoir.mass_properties.mass
        weight.thermal_management_system.total = tms_weight

        # -------------------------------------------------------------------------------
        # Wing and empennage — CADDEE_LPC nasa_lpc regression (replaces EVTOL.compute_wing_weight)
        # -------------------------------------------------------------------------------
        maxSpan = 0
        for wing in vehicle.wings:
            maxSpan = max(wing.spans.projected, maxSpan)

            if wing.symbolic:
                wing_weight = 0
                weight.wings[wing.tag] = wing_weight
                wing.mass_properties.mass = wing_weight
                weight.wings_total += wing_weight
            elif isinstance(wing, RCAIDE.Library.Components.Wings.Main_Wing):
                wing_weight = CADDEE_LPC.compute_wing_weight(
                    wing_area=wing_area_for_regression, wing_AR=wing_AR_for_regression,
                    fuselage_length=fuselage_length_for_regression,
                    battery_mass=battery_mass_design_var, cruise_speed=cruise_speed,
                    correction_factor=correction_factor,
                )
                weight.wings[wing.tag] = wing_weight
                wing.mass_properties.mass = wing_weight
                weight.wings_total += wing_weight
            elif isinstance(wing, (RCAIDE.Library.Components.Wings.Horizontal_Tail,
                                    RCAIDE.Library.Components.Wings.Vertical_Tail)):
                # handled below as a single combined empennage regression, not per-surface
                pass
            else:
                wing_weight = EVTOL.compute_wing_weight(wing, vehicle, maxLift/5, safety_factor=safety_factor, max_g_load=max_g_load)
                weight.wings[wing.tag] = wing_weight
                wing.mass_properties.mass = wing_weight
                weight.wings_total += wing_weight

            if isinstance(wing, RCAIDE.Library.Components.Wings.Main_Wing):
                wiring_weight = EVTOL_Common.compute_wiring_weight(wing, vehicle, maxLiftPower/(eta*total_number_of_rotors)) * Units.kg
            else:
                wiring_weight = 0
            weight.wiring += wiring_weight

        h_tail_area = 0.0
        v_tail_area = 0.0
        for wing in vehicle.wings:
            if isinstance(wing, RCAIDE.Library.Components.Wings.Horizontal_Tail):
                h_tail_area += wing.areas.reference
            elif isinstance(wing, RCAIDE.Library.Components.Wings.Vertical_Tail):
                v_tail_area += wing.areas.reference

        if h_tail_area > 0 or v_tail_area > 0:
            empennage_weight = CADDEE_LPC.compute_empennage_weight(
                h_tail_area=h_tail_area, v_tail_area=v_tail_area, correction_factor=correction_factor,
            )
            for wing in vehicle.wings:
                if isinstance(wing, (RCAIDE.Library.Components.Wings.Horizontal_Tail,
                                      RCAIDE.Library.Components.Wings.Vertical_Tail)):
                    share = (wing.areas.reference / (h_tail_area + v_tail_area)) if (h_tail_area + v_tail_area) > 0 else 0.0
                    wing.mass_properties.mass = empennage_weight * share
            weight.empennage_total = empennage_weight

        # -------------------------------------------------------------------------------
        # Landing Gear Weight — untouched from Physics_Based
        # -------------------------------------------------------------------------------
        for LG in vehicle.landing_gears:
            if isinstance(LG, RCAIDE.Library.Components.Landing_Gear.Main_Landing_Gear):
                LG.mass_properties.mass = 2 / 3 * weight.landing_gear
            if isinstance(LG, RCAIDE.Library.Components.Landing_Gear.Nose_Landing_Gear):
                LG.mass_properties.mass = 1 / 3 * weight.landing_gear

        # -------------------------------------------------------------------------------
        # Fuselage Weight — CADDEE_LPC nasa_lpc regression (replaces EVTOL.compute_fuselage_weight)
        # -------------------------------------------------------------------------------
        for fuse in vehicle.fuselages:
            fuselage_weight = CADDEE_LPC.compute_fuselage_weight(
                wing_area=wing_area_for_regression, wing_AR=wing_AR_for_regression,
                fuselage_length=fuselage_length_for_regression,
                battery_mass=battery_mass_design_var, cruise_speed=cruise_speed,
                correction_factor=correction_factor,
            )
            fuse.mass_properties.center_of_gravity[0][0] = .45*fuse.lengths.total
            fuse.mass_properties.mass = fuselage_weight
            weight.fuselage += fuselage_weight

        # -------------------------------------------------------------------------------
        # Boom Weight — CADDEE_LPC nasa_lpc regression (replaces EVTOL.compute_boom_weight)
        # -------------------------------------------------------------------------------
        n_booms = len(list(vehicle.booms))
        if n_booms:
            total_boom_weight = CADDEE_LPC.compute_boom_weight(
                wing_area=wing_area_for_regression, wing_AR=wing_AR_for_regression,
                fuselage_length=fuselage_length_for_regression,
                battery_mass=battery_mass_design_var, cruise_speed=cruise_speed,
                correction_factor=correction_factor,
            )
            for boom in vehicle.booms:
                boom_weight = total_boom_weight / n_booms
                weight.booms += boom_weight
                boom.mass_properties.mass = boom_weight

        # -------------------------------------------------------------------------------
        # Pack Up Outputs — same schema as Physics_Based
        # -------------------------------------------------------------------------------
        output = Data()
        output.empty = Data()

        output.empty.structural = Data()
        output.empty.structural.wings = miscelleneous_weight_factor * weight.wings_total
        output.empty.structural.fuselage = miscelleneous_weight_factor * weight.fuselage
        output.empty.structural.empennage = miscelleneous_weight_factor * weight.empennage_total
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
        output.empty.systems.control_systems = 0.0
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
