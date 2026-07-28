# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/Chakraborty/compute_operating_empty_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

import RCAIDE
from RCAIDE.Framework.Core import Units, Data
import RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.Physics_Based as Vahana
import RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.Common as EVTOL_Common
import RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Conventional.General_Aviation.Raymer as Raymer
import RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Conventional.General_Aviation.Roskam as Roskam
from RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Conventional.Common import compute_payload_weight

import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  Compute Operating Empty Weight — Chakraborty method
# ----------------------------------------------------------------------------------------------------------------------
# Direct copy of `Physics_Based/compute_operating_empty_weight.py` — same pattern already used for
# `NDARC`/`Hydra`/`CADDEE_LPC` (a copy with component calls swapped in-place, inside the same
# converging loop, not a separate decoupled pass — see the inventory doc's 2026-07-28 lesson on
# why the decoupled-pass approach is wrong).
#
# Ports Chakraborty & Mishra, "Generalized Energy-Based Flight Vehicle Sizing and Performance
# Analysis Methodology," Journal of Aircraft, Vol. 58, No. 4, 2021, Appendix C's airframe-
# structure methodology: General Aviation statistical weight-estimation relationships (WERs) from
# Raymer and Roskam, averaged per component where the paper cites both, with a tilt-mechanism
# weight penalty (`KTW`, default 1.20) applied to wing/horizontal-tail/canard. See
# `20-rcaide-weight-method-porting-inventory.md` for the full Table C1 breakdown, the Raymer-
# textbook cross-check, and the 3 real RCAIDE bugs found and fixed along the way (unrelated to
# this port itself).
#
# **Not physics-based** — these are statistical regressions fit to existing GA aircraft, not a
# structural sizing method. Included specifically because the user wanted an independent
# statistical cross-check for the standing "why does Vahana run light" question, not because this
# method is expected to be more physically correct than `Vahana`/`NDARC`/`Hydra`.
#
# Per-component method map (source: paper's own Table C1 attribution):
#   Fuselage           avg(Roskam Wfus,1, Raymer Wfus,2)
#   Wing               avg(Roskam Wwing,1, Roskam Wwing,2, Raymer Wwing,3), each x KTW
#   Horizontal tail     Raymer Wht x KTW
#   Vertical tail       avg(Roskam Wvt,1, Raymer Wvt,2)
#   Landing gear        Raymer only (main + nose), no averaging
#   Electrical          avg(Roskam Welec,1, Raymer Welec,2)
#   Fuel system         Raymer only (via Raymer's bundled compute_systems_weight)
#   Flight controls     Raymer only (via the same bundled call)
#   Furnishings         avg(Roskam Wfur,1, Roskam Wfur,2) -- paper cites both furnishings variants
#                       to Roskam; no Raymer furnishings WER exists in Table C1 at all
#   Hydraulics          Raymer only (via the same bundled call, now source-corrected)
#
# Components with no Raymer/Roskam WER in the paper's table at all (motor, rotor/blade, BRS,
# battery, thermal management, wiring, seats, avionics-as-a-component, ECS, propeller [the
# paper's own propeller source is Plencner, NASA TM-83458 -- not ported]) fall back to `Vahana`,
# same role it plays for `NDARC`/`Hydra`/`CADDEE_LPC`'s own gaps.
#
# **New per-vehicle input convention introduced here** — a `.Chakraborty` namespace:
#   vehicle.Chakraborty.KTW              float, default 1.20 (the paper's own tilt-mechanism
#                                         penalty). Set to 1.0 for a conventional (non-tilting)
#                                         fixed-wing GA aircraft.
#   vehicle.Chakraborty.fuel_volume      float, default 0.0 m^3 -- total fuel volume (Qtot in the
#                                         paper's own notation). 0.0 for an all-electric vehicle.
#   vehicle.Chakraborty.internal_fuel_volume
#                                         float, default = fuel_volume -- internal/integral fuel
#                                         volume (Qint). Only matters when some fuel is carried in
#                                         non-integral (e.g. external/droppable) tanks; defaulting
#                                         to all-internal is the common case.
#   vehicle.Chakraborty.number_of_fuel_tanks / number_of_engines
#                                         int, default 1 each -- only meaningful if fuel_volume>0.


def compute_operating_empty_weight(vehicle, settings=None):
    """ Calculates the empty vehicle mass for a (typically tiltwing) eVTOL aircraft using General
        Aviation statistical weight-estimation relationships from Raymer and Roskam (averaged per
        component, per Chakraborty & Mishra 2021's own methodology), for airframe structure, and
        `Vahana`/`Physics_Based`'s own methods for everything else (motor, rotor, BRS, battery,
        seats/avionics/ECS, wiring, thermal management) -- neither Raymer nor Roskam's GA WERs
        cover any of those. See module docstring for the full method map, the `.Chakraborty`
        namespace convention, and the "not physics-based" caveat.

        Inputs:
            vehicle:     RCAIDE Vehicle Data Structure (with `.Chakraborty` namespace inputs set
                         -- see module docstring)
            settings:    Weights analysis settings (same object every other method uses;
                         `settings.miscelleneous_weight_factor` applied identically)

        Outputs:
            output:      Data dictionary, same schema as `Physics_Based.compute_operating_empty_weight`
    """
    v_chak = getattr(vehicle, 'Chakraborty', Data())
    KTW = getattr(v_chak, 'KTW', 1.20)
    fuel_volume = getattr(v_chak, 'fuel_volume', 0.0)
    internal_fuel_volume = getattr(v_chak, 'internal_fuel_volume', fuel_volume)
    number_of_fuel_tanks = getattr(v_chak, 'number_of_fuel_tanks', 1)
    number_of_engines = getattr(v_chak, 'number_of_engines', 1)

    main_wing = None
    for wing in vehicle.wings:
        if isinstance(wing, RCAIDE.Library.Components.Wings.Main_Wing):
            main_wing = wing
    if main_wing is None:
        raise ValueError(
            "Chakraborty.compute_operating_empty_weight: no Main_Wing found on vehicle."
        )

    ref_fuselage = None
    for fuse in vehicle.fuselages:
        ref_fuselage = fuse
        break
    if ref_fuselage is None:
        raise ValueError(
            "Chakraborty.compute_operating_empty_weight: no fuselage found on vehicle."
        )

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
        weight.ECS = vehicle.number_of_passengers * 7. * Units.kg

        # -------------------------------------------------------------------------------
        # Network Weight — rotor/motor/BRS, untouched from Physics_Based (no Chakraborty WER
        # covers propulsion at all)
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
                    maxVTip = rotor.hover.design_angular_velocity * rotor.tip_radius
                    lift_rotor_servo_weight = 0.65 * Units.kg
                    if rotor.oei.design_thrust is None:
                        design_thrust = rotor.hover.design_thrust
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
        # Raymer's bundled systems call — gives real fuel-system/flight-controls/hydraulics
        # (all now source-corrected), plus intermediate avionics/fuel-system values the
        # Chakraborty electrical average needs. Furnishings/electrical from this call are
        # discarded and replaced below with the Roskam-averaged versions.
        # -------------------------------------------------------------------------------
        raymer_systems = Raymer.compute_systems_weight(
            vehicle, V_fuel=fuel_volume, V_int=internal_fuel_volume,
            N_tank=number_of_fuel_tanks, N_eng=number_of_engines,
        )

        # -------------------------------------------------------------------------------
        # Wing, horizontal tail — averaged Raymer+Roskam, KTW applied
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
                w_raymer = Raymer.compute_main_wing_weight(wing, vehicle, m_fuel=0.0) * KTW
                w_roskam_1 = Roskam.compute_wing_weight_1(vehicle, wing, KTW=KTW)
                w_roskam_2 = Roskam.compute_wing_weight_2(vehicle, wing, KTW=KTW)
                wing_weight = (w_raymer + w_roskam_1 + w_roskam_2) / 3.0
                weight.wings[wing.tag] = wing_weight
                wing.mass_properties.mass = wing_weight
                weight.wings_total += wing_weight
            elif isinstance(wing, RCAIDE.Library.Components.Wings.Horizontal_Tail):
                wing_weight = Raymer.compute_horizontal_tail_weight(wing, vehicle) * KTW
                weight.wings[wing.tag] = wing_weight
                wing.mass_properties.mass = wing_weight
                weight.wings_total += wing_weight
            elif isinstance(wing, RCAIDE.Library.Components.Wings.Vertical_Tail):
                w_raymer_vt = Raymer.compute_vertical_tail_weight(wing, vehicle)
                w_roskam_vt = Roskam.compute_vertical_tail_weight(wing, vehicle)
                wing_weight = (w_raymer_vt + w_roskam_vt) / 2.0
                weight.wings[wing.tag] = wing_weight
                wing.mass_properties.mass = wing_weight
                weight.wings_total += wing_weight
            else:
                wing_weight = Vahana.compute_wing_weight(wing, vehicle, maxLift/5, safety_factor=safety_factor, max_g_load=max_g_load)
                weight.wings[wing.tag] = wing_weight
                wing.mass_properties.mass = wing_weight
                weight.wings_total += wing_weight

            if isinstance(wing, RCAIDE.Library.Components.Wings.Main_Wing):
                wiring_weight = EVTOL_Common.compute_wiring_weight(wing, vehicle, maxLiftPower/(eta*total_number_of_rotors) if eta and total_number_of_rotors else 0.0) * Units.kg
            else:
                wiring_weight = 0
            weight.wiring += wiring_weight

        # -------------------------------------------------------------------------------
        # Landing Gear — Raymer only, real Wmlg/Wnlg formulas (not Vahana's fixed 2% fraction)
        # -------------------------------------------------------------------------------
        strut_length_main = 0.0
        strut_length_nose = 0.0
        for LG in vehicle.landing_gears:
            if isinstance(LG, RCAIDE.Library.Components.Landing_Gear.Main_Landing_Gear):
                strut_length_main = getattr(LG, 'strut_length', 0.3 * Units.ft)
            if isinstance(LG, RCAIDE.Library.Components.Landing_Gear.Nose_Landing_Gear):
                strut_length_nose = getattr(LG, 'strut_length', 0.2 * Units.ft)
        if strut_length_main == 0.0:
            strut_length_main = 0.3 * Units.ft
        if strut_length_nose == 0.0:
            strut_length_nose = 0.2 * Units.ft

        lg = Raymer.compute_landing_gear_weight(MTOW, Nult=vehicle.flight_envelope.ultimate_load,
                                                 strut_length_main=strut_length_main,
                                                 strut_length_nose=strut_length_nose)
        weight.landing_gear = lg.main + lg.nose
        for LG in vehicle.landing_gears:
            if isinstance(LG, RCAIDE.Library.Components.Landing_Gear.Main_Landing_Gear):
                LG.mass_properties.mass = lg.main
            if isinstance(LG, RCAIDE.Library.Components.Landing_Gear.Nose_Landing_Gear):
                LG.mass_properties.mass = lg.nose

        # -------------------------------------------------------------------------------
        # Fuselage — averaged Raymer+Roskam
        # -------------------------------------------------------------------------------
        l_ht = 0.4 * ref_fuselage.lengths.total
        for wing in vehicle.wings:
            if isinstance(wing, RCAIDE.Library.Components.Wings.Horizontal_Tail):
                l_ht = abs(wing.origin[0][0] - main_wing.origin[0][0])

        f_raymer = Raymer.compute_fuselage_weight(ref_fuselage, vehicle, l_ht=l_ht)
        f_roskam = Roskam.compute_fuselage_weight(ref_fuselage, vehicle)
        fuselage_weight = (f_raymer + f_roskam) / 2.0
        ref_fuselage.mass_properties.center_of_gravity[0][0] = .45 * ref_fuselage.lengths.total
        ref_fuselage.mass_properties.mass = fuselage_weight
        weight.fuselage = fuselage_weight

        # -------------------------------------------------------------------------------
        # Electrical, furnishings — averaged Raymer+Roskam where the paper does; fuel system/
        # flight controls/hydraulics straight from Raymer's bundled call above
        # -------------------------------------------------------------------------------
        elec_roskam = Roskam.compute_electrical_weight(raymer_systems.W_fuel_system, raymer_systems.W_avionics)
        elec_raymer2 = Raymer.compute_electrical_weight(vehicle)
        electrical_weight = (elec_roskam + elec_raymer2) / 2.0

        furn_1 = Roskam.compute_furnishings_weight_1(vehicle)
        furn_2 = Roskam.compute_furnishings_weight_2(vehicle)
        furnishings_weight = (furn_1 + furn_2) / 2.0

        # -------------------------------------------------------------------------------
        # Boom (no Chakraborty/Raymer/Roskam GA model — Vahana/RCAIDE only)
        # -------------------------------------------------------------------------------
        for boom in vehicle.booms:
            boom_weight = Vahana.compute_boom_weight(boom) * Units.kg
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
        output.empty.propulsion.fuel_system = miscelleneous_weight_factor * raymer_systems.W_fuel_system
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
        output.empty.systems.control_systems = miscelleneous_weight_factor * raymer_systems.W_flight_control
        output.empty.systems.apu = 0.0
        output.empty.systems.electrical = miscelleneous_weight_factor * electrical_weight
        output.empty.systems.avionics = miscelleneous_weight_factor * weight.avionics
        output.empty.systems.hydraulics = miscelleneous_weight_factor * raymer_systems.W_hyd_pnu
        output.empty.systems.furnishings = miscelleneous_weight_factor * furnishings_weight
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
