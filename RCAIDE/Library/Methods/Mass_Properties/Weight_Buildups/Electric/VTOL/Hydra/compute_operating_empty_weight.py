# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/Hydra/compute_operating_empty_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

# RCAIDE
import RCAIDE
from RCAIDE.Framework.Core import Units, Data
import RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.Hydra as Hydra
import RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.NDARC as NDARC
import RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.Physics_Based as Vahana
import RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.Common as EVTOL_Common
from RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Conventional.Common import compute_payload_weight

# package imports
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  Compute Operating Empty Weight — Hydra method
# ----------------------------------------------------------------------------------------------------------------------
# Orchestrates every validated Hydra component (see
# 01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md for the
# component-by-component porting/validation record) into the same output schema
# `Physics_Based/NDARC` return, so downstream RCAIDE machinery (CG, MOI, AVL mass file, plots)
# works identically regardless of which method is selected.
#
# **New per-vehicle input convention introduced here** — a `.Hydra` namespace on wings/rotors,
# mirroring `.NDARC`'s own namespace. None of the three current test vehicles set these yet —
# wiring a real vehicle through `method='Hydra'` is a separate follow-up step, exactly the same
# precedent as `NDARC`'s orchestrator. Required/optional fields, by object:
#
#   rotor.Hydra.material         str, REQUIRED — spar material name (see
#                                 `Hydra/compute_rotor_weight.py::_spar_properties` for valid
#                                 values), OR a real `Solid` material instance (e.g.
#                                 `AS4_3502_Unidirectional_Carbon_Fiber`), 2026-07-30 — only
#                                 `.density`/`.ultimate_tensile_strength` are consumed either way.
#                                 No default assumed — blade mass is highly sensitive to
#                                 material choice, same reasoning as `NDARC`'s `nu_blade`.
#   rotor.Hydra.load_factor      float, REQUIRED — limit load factor (nz). No default assumed,
#                                 same reasoning as material — this drives the entire spanwise spar
#                                 sizing loop.
#   rotor.Hydra.precone_deg      float, default 3.0 (matches source default).
#   rotor.Hydra.rho_filler       float, default 52.0 kg/m^3 (matches source default).
#   rotor.Hydra.tech_factor      float, default 1.0 — uniform scale on blade+hub+actuator mass.
#                                 Exposed 2026-07-29 after calibrating against two real EASA Type
#                                 Certificate Data Sheets (Hoffmann HO-V62, E-Props EPGU3/
#                                 GLORIEUSE-3) implied roughly 0.2-0.8 for blade+hub specifically,
#                                 while the actuator term alone checked out close to 1.0 already —
#                                 see 01-mission-profiles/00-decisions/2026-07-29-rotor-motor-
#                                 weight-formula-comparison.md. Left at the default 1.0 here
#                                 (unvalidated end-to-end); a single scalar cannot correct
#                                 blade+hub and actuator independently, so treat any non-default
#                                 value as a rough interim correction, not a final calibration.
#
#   wing.Hydra.tip_propulsor_tags   list[str], same convention as `NDARC`'s field of the same
#                                    name — which propulsors belong to this wing's *tilting* group
#                                    (whirl-flutter-relevant, per NDARC's own Theory Manual
#                                    convention). Do NOT use this for boom-mounted, non-tilting
#                                    rotors that still load the wing structurally — see
#                                    `load_bearing_stations` below, a deliberately separate field
#                                    since the two questions ("does this rotor tilt" vs "does this
#                                    rotor load the wing") are physically different (2026-07-28).
#   wing.Hydra.load_bearing_stations   list[dict], each {'y': spanwise position [m],
#                                    'propulsor_tags': [str, ...], 'boom_tag': str or None} — one
#                                    entry per lumped mass station the wing structurally carries
#                                    (e.g. one entry per boom for a boom-mounted-rotor vehicle like
#                                    lift_cruise/tilt_stopped_rotor_v_tail). Real converged
#                                    rotor+motor mass for the tagged propulsors, plus the tagged
#                                    boom's own structural mass (if `boom_tag` given), are summed
#                                    into that station's lumped mass. Takes priority over
#                                    `tip_propulsor_tags` if both are set (should not normally both
#                                    be set on the same wing).
#   wing.Hydra.lift_fraction        float, default 1.0 — fraction of MTOW-derived lift this wing
#                                    group carries (source: `group.lift_frac/group.nwings`).
#   wing.Hydra.spar_material        `Solid` instance, optional -- 2026-07-30. Drives the spar/
#                                    motor-mount sizing. Unset defaults to a real, sourced
#                                    AS4_3502_Unidirectional_Carbon_Fiber instance (UNNOTCHED
#                                    B-basis -- the wing spar is a continuous laminate beam, no
#                                    fastener modeled, so open-hole data doesn't apply here;
#                                    `compute_wing_weight.py::_default_spar_material()`) -- NOT the
#                                    original hardcoded 122 GPa/1650 kg/m^3/275 MPa/47 MPa
#                                    constants (a real default-behavior change, not backward-
#                                    compatible; see the 2026-07-30 decision doc). Only used on
#                                    load-bearing-station/tip-tag (rotor-carrying) wings -- the
#                                    Vahana fallback branch (non-rotor-carrying wings/tails) is
#                                    unaffected, out of scope this session.
#
#   fuselage.keel_materials.root_bending_moment_carrier   `Solid` instance, optional -- 2026-07-30.
#                                    Same field `Vahana`/`Physics_Based` already reads (not a new
#                                    `.Hydra`-specific namespace) -- drives the keel bending term's
#                                    allowable stress/density in both methods identically. Unset
#                                    falls back to Hydra's own hardcoded 450 MPa/1660 kg/m^3.
#   fuselage.keel_materials.shear_carrier         `Solid` instance, optional -- same convention,
#                                    drives keel torsion/shear term. Unset falls back to 47 MPa.
#   fuselage.keel_materials.bearing_carrier       `Solid` instance, optional -- drives keel
#                                    landing-bearing term. Unset falls back to 400 MPa.
#   fuselage.materials.bolt_materials.landing_pad_bolt   `Solid` instance, optional -- drives bolt
#                                    shear sizing. Unset falls back to 500 MPa.
#
# Components with no Hydra-specific model (fall back to Vahana/Common, matching the porting
# inventory's own component matrix): flight control system (confirmed dead in source, see
# inventory doc — always 0.0), boom, wiring, payload, battery/fuel-cell pass-through, thermal
# management, and seats/avionics/ECS (same flat-constant fallback as `Vahana`/`NDARC`, applies
# identically here — see the inventory doc's "Avionics / seats / ECS / systems" row).


def compute_operating_empty_weight(vehicle, settings=None):
    """ Calculates the empty vehicle mass for an eVTOL-type aircraft using Hydra-derived component
        weight methods where a validated Hydra port exists, falling back to `Vahana`/`Physics_Based`
        methods elsewhere. See module docstring for the full per-component method map and the new
        `.Hydra` namespace input convention this orchestrator introduces.

        Inputs:
            vehicle:     RCAIDE Vehicle Data Structure (with `.Hydra` namespace inputs set — see
                         module docstring)
            settings:    Weights analysis settings (same object `Physics_Based`/`NDARC` use;
                         `settings.miscelleneous_weight_factor` applied identically)

        Outputs:
            output:      Data dictionary, same schema as `Physics_Based.compute_operating_empty_weight`
    """
    diff = 100
    iterations = 0
    MTOW = vehicle.mass_properties.max_takeoff
    tolerance = settings.mtow_convergence_tolerance

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
        # Flight controls, 2026-08-10: was hardcoded 0.0. CORRECTION (same day, after this was
        # first wired in): the 0.0 was actually a verified-correct fidelity match to Hydra's own
        # real source, not a gap - the porting inventory doc (01-docs/rcaide/20-...md, based on a
        # direct 2026-07-28 source read) confirms Hydra's `afdd/flight_controls.py::flightctrl_
        # weight` is imported but never called anywhere in Hydra's own code; the real
        # `vehicle_empty_weight()` return dict has no flight-control key at all. So
        # `WeightSettings.flight_controls_constant=11.5`/`_exponent=0.4` (comment: "from hydra")
        # is MISLABELED - it is not Hydra's formula (Hydra has none); its real origin looks like
        # a separate generic GA-style regression already used in AircraftDesign's own legacy
        # `uav_design.py`/`impulls.py` tools. Deliberately wired in here anyway (user decision,
        # 2026-08-10) as a practical non-zero estimate, accepting the departure from pure Hydra
        # fidelity for this one component - zero is obviously wrong for a real vehicle. Cross-
        # checked against NDARC's own simplest fixed-wing flight-controls fraction model
        # (0.91*MTOW_lbs^0.6, same MTOW-power-law structure) for these 3 vehicles' real MTOWs -
        # NDARC's analog came out roughly 7x higher, so this is likely OPTIMISTIC (under-
        # predicts) - not a final calibration. See 01-mission-profiles/00-decisions/
        # 2026-08-10-o40-registry-scope-expanded-full-audit-needed.md.
        #
        # `settings` here is RCAIDE's own Electric_VTOL().settings (a bare Data() object,
        # default fields just mtow_convergence_tolerance/miscelleneous_weight_factor) - NOT
        # AircraftDesign's WeightSettingsRCAIDEHydra, confirmed 2026-08-10 by directly running
        # the real dispatch path (debug_hydra_orchestrator_test.py), which crashed with
        # AttributeError on a first cut of this code that assumed the settings values were
        # already live-wired through. Every field below uses getattr(...) with the
        # WeightSettingsRCAIDEHydra-confirmed default as a fallback, same idiom this file
        # already uses for rotor.Hydra.tech_factor etc. - picks up a real WeightSettingsRCAIDEHydra
        # value if/when one is ever actually attached to weights_analysis.settings (not done
        # anywhere today - a real, deeper gap than the registry's "already wired" notes assumed),
        # falls back to the confirmed value otherwise.
        flight_controls_constant = getattr(settings, 'flight_controls_constant', 11.5)
        flight_controls_exponent = getattr(settings, 'flight_controls_exponent', 0.4)
        mass_in_lbs = MTOW / Units.lbs
        weight.flight_controls = (flight_controls_constant
                                   * (mass_in_lbs / 1000.) ** flight_controls_exponent
                                   * Units.lbs)
        weight.thermal_management_system = Data()

        # -------------------------------------------------------------------------------
        # Payload
        # -------------------------------------------------------------------------------
        # has_pilot, 2026-08-10: pilot is tracked separately from vehicle.number_of_passengers
        # (confirmed by user 2026-08-10: number_of_passengers is passenger-only, does not
        # already include a pilot) - adds one pilot's own weight (as flight crew, below) plus
        # one extra seat and ECS share when a vehicle carries one. No default assumed (getattr
        # falls back to False) - a vehicle must set this explicitly, same reasoning as
        # rotor.Hydra.material/load_factor having no silent default.
        n_pilots = 1 if getattr(vehicle, 'has_pilot', False) else 0
        passenger_weight = getattr(settings, 'passenger_weight', 100.)
        seat_weight = getattr(settings, 'seat_weight', 15.)
        avionics_mass_kg = getattr(settings, 'avionics_mass_kg', 15.)
        ecs_per_passenger = getattr(settings, 'ecs_per_passenger', 7.)
        payload = compute_payload_weight(vehicle, W_passenger=passenger_weight * Units.kg,
                                          W_baggage=0 * Units.lbs)
        weight.seats = (vehicle.number_of_passengers + n_pilots) * seat_weight * Units.kg
        weight.passengers = payload.passengers
        weight.avionics = avionics_mass_kg * Units.kg
        weight.landing_gear = 0.0
        weight.ECS = (vehicle.number_of_passengers + n_pilots) * ecs_per_passenger * Units.kg

        # -------------------------------------------------------------------------------
        # Rotors + motors: blade/hub/actuator and motor via Hydra's live DC_motor/blade_wt_modelv2
        # -------------------------------------------------------------------------------
        total_number_of_rotors = 0
        propulsor_data = {}
        for network in vehicle.networks:
            for propulsor in network.propulsors:
                rotor = propulsor.rotor
                r_hydra = getattr(rotor, 'Hydra', Data())
                material = getattr(r_hydra, 'material', None)
                load_factor = getattr(r_hydra, 'load_factor', None)
                if material is None:
                    raise ValueError(
                        f"Hydra.compute_operating_empty_weight: rotor '{rotor.tag}' has no "
                        "rotor.Hydra.material set (spar material, required — no default assumed, "
                        "see module docstring)."
                    )
                if load_factor is None:
                    raise ValueError(
                        f"Hydra.compute_operating_empty_weight: rotor '{rotor.tag}' has no "
                        "rotor.Hydra.load_factor set (limit load factor, required — no default "
                        "assumed, see module docstring)."
                    )
                precone_deg = getattr(r_hydra, 'precone_deg', 3.0)
                rho_filler = getattr(r_hydra, 'rho_filler', 52.0)
                tech_factor = getattr(r_hydra, 'tech_factor', 0.5)

                n_blade = rotor.number_of_blades
                radius = rotor.tip_radius
                chord = np.mean(rotor.chord_distribution)
                omega = EVTOL_Common.max_design_field(rotor, 'design_angular_velocity')
                thrust = EVTOL_Common.max_design_field(rotor, 'design_thrust')

                blade_hub, rotor_total = Hydra.compute_rotor_weight(
                    radius=radius, chord=chord, omega=omega, thrust=thrust, n_blade=n_blade,
                    n_rotor=1, material=material, load_factor=load_factor,
                    precone_deg=precone_deg, rho_filler=rho_filler, tech_factor=tech_factor,
                )
                weight.rotors += blade_hub['blades']
                weight.hubs += blade_hub['hub'] + blade_hub['actuator']
                rotor.mass_properties.mass = rotor_total

                torque = propulsor.motor.design_torque
                # 2026-07-30: ported from Hydra's own DC_motor formula to NDARC's empirical_2026
                # regression (real 30-motor fit, mean|err|=19.0%) -- motor method decided
                # independently of the fuselage/wing/rotor chassis choice, see
                # 00-decisions/2026-07-29-rotor-motor-weight-formula-comparison.md. `Hydra`'s own
                # `compute_motor_weight` (DC_motor-based) is no longer called here.
                motor_mass = NDARC.compute_motor_weight(design_torque=torque, method='empirical_2026')
                weight.motors += motor_mass
                propulsor.motor.mass_properties.mass = motor_mass

                total_number_of_rotors += 1
                propulsor_data[propulsor.tag] = {
                    'y': abs(rotor.origin[0][1]),
                    'rotor_mass_no_mount': motor_mass + rotor_total,
                    'thrust': thrust,
                    'radius': radius,
                }

        weight.BRS = Hydra.compute_emergency_system_weight(MTOW, tech_factor=1.0)

        lg = Hydra.compute_landing_gear_weight(MTOW, tech_factor=1.0)
        weight.landing_gear = lg['total']

        # -------------------------------------------------------------------------------
        # Thermal management (pass-through, same as Vahana/NDARC — no Hydra model)
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
        # Battery / fuel cell (pass-through, same as Vahana/NDARC)
        # -------------------------------------------------------------------------------
        for network in vehicle.networks:
            for bus in network.busses:
                for modules in bus.battery_modules:
                    weight.battery += modules.mass_properties.mass * Units.kg
                for fuel_cell in bus.fuel_cell_stacks:
                    weight.fuel_cell += fuel_cell.mass_properties.mass * Units.kg

        # -------------------------------------------------------------------------------
        # Wings and tails — Hydra target-frequency spar model if rotor-carrying, Vahana
        # beam-bending otherwise
        # -------------------------------------------------------------------------------
        maxLift = MTOW * 1.1 * 9.81
        maxSpan = 0
        wing_lift_fractions = []
        for wing in vehicle.wings:
            maxSpan = max(wing.spans.projected, maxSpan)
            w_hydra = getattr(wing, 'Hydra', Data())
            tip_tags = getattr(w_hydra, 'tip_propulsor_tags', [])
            load_bearing_stations = getattr(w_hydra, 'load_bearing_stations', [])
            lift_fraction = getattr(w_hydra, 'lift_fraction', 1.0)
            spar_material = getattr(w_hydra, 'spar_material', None)
            wing_lift_fractions.append(lift_fraction)

            if wing.symbolic:
                wing_weight = 0
            elif load_bearing_stations:
                y_positions = [station['y'] for station in load_bearing_stations]
                rotor_masses_arr = []
                rotor_thrusts_arr = []
                for station in load_bearing_stations:
                    station_mass = sum(propulsor_data[tag]['rotor_mass_no_mount'] for tag in station['propulsor_tags'])
                    station_thrust = sum(propulsor_data[tag]['thrust'] for tag in station['propulsor_tags'])
                    boom_tag = station.get('boom_tag')
                    if boom_tag is not None:
                        boom = next(b for b in vehicle.booms if b.tag == boom_tag)
                        station_mass += Vahana.compute_boom_weight(boom) * Units.kg
                    rotor_masses_arr.append(station_mass)
                    rotor_thrusts_arr.append(station_thrust)
                ref_tag = load_bearing_stations[0]['propulsor_tags'][0]
                ref = propulsor_data[ref_tag]

                group = Hydra.compute_wing_weight_group(
                    vehicle_mtow=MTOW, n_wings=1, aspect_ratio=wing.aspect_ratio,
                    area=wing.areas.reference, lift_fraction=lift_fraction,
                    rotor_radius=ref['radius'], rotor_mass_assembly=ref['rotor_mass_no_mount'],
                    rotor_y_positions=y_positions, rotor_masses=rotor_masses_arr,
                    rotor_thrusts=rotor_thrusts_arr, tech_factor_wing=1.0,
                    tech_factor_flight_control=1.0, spar_material=spar_material,
                )
                wing_weight = group['structure'] + group['actuators'] + group['tilters'] + group['mounts']
            elif tip_tags:
                y_positions = [propulsor_data[tag]['y'] for tag in tip_tags]
                rotor_masses_arr = [propulsor_data[tag]['rotor_mass_no_mount'] for tag in tip_tags]
                rotor_thrusts_arr = [propulsor_data[tag]['thrust'] for tag in tip_tags]
                ref = propulsor_data[tip_tags[0]]

                group = Hydra.compute_wing_weight_group(
                    vehicle_mtow=MTOW, n_wings=1, aspect_ratio=wing.aspect_ratio,
                    area=wing.areas.reference, lift_fraction=lift_fraction,
                    rotor_radius=ref['radius'], rotor_mass_assembly=ref['rotor_mass_no_mount'],
                    rotor_y_positions=y_positions, rotor_masses=rotor_masses_arr,
                    rotor_thrusts=rotor_thrusts_arr, tech_factor_wing=1.0,
                    tech_factor_flight_control=1.0, spar_material=spar_material,
                )
                wing_weight = group['structure'] + group['actuators'] + group['tilters'] + group['mounts']
            else:
                wing_weight = Vahana.compute_wing_weight(wing, vehicle, maxLift / 5)

            wing_tag = wing.tag
            weight.wings[wing_tag] = wing_weight
            wing.mass_properties.mass = wing_weight
            weight.wings_total += wing_weight

            if isinstance(wing, RCAIDE.Library.Components.Wings.Main_Wing):
                # cablePower left at 0.0 — matches NDARC's own orchestrator exactly: its
                # `maxLiftPower/(eta*total_number_of_rotors)` expression evaluates to 0.0 too,
                # since `maxLiftPower` is initialized but never updated in that file either. Not
                # fixed here — a real gap in both orchestrators, flagged in the inventory doc as a
                # shared follow-up, out of scope for finishing this method.
                wiring_weight = EVTOL_Common.compute_wiring_weight(
                    wing, vehicle, 0.0
                ) * Units.kg
            else:
                wiring_weight = 0
            weight.wiring += wiring_weight

        min_wing_lift_fraction = min(wing_lift_fractions) if wing_lift_fractions else 1.0

        # -------------------------------------------------------------------------------
        # Fuselage — Hydra "GB/ZL" ellipsoid model
        # -------------------------------------------------------------------------------
        for fuse in vehicle.fuselages:
            # 2026-07-30: reuse the SAME `fuse.keel_materials`/`fuse.materials.bolt_materials`
            # real-material-object convention `Vahana` already reads (Physics_Based/
            # compute_fuselage_weight.py) -- one vehicle-level material input now drives both
            # methods identically, rather than each method needing its own namespace. Any field
            # left unset falls back to this module's original hardcoded constants (see
            # `compute_fuselage_weight`'s `keel_materials` docstring) -- bit-for-bit unchanged
            # when `fuse.keel_materials` is not set, exactly like every other vehicle here today.
            keel_mats = getattr(fuse, 'keel_materials', None)
            bolt_mats = getattr(getattr(fuse, 'materials', None), 'bolt_materials', None)
            rbm_mat = getattr(keel_mats, 'root_bending_moment_carrier', None) if keel_mats else None
            shear_mat = getattr(keel_mats, 'shear_carrier', None) if keel_mats else None
            bearing_mat = getattr(keel_mats, 'bearing_carrier', None) if keel_mats else None
            bolt_mat = getattr(bolt_mats, 'landing_pad_bolt', None) if bolt_mats else None

            hydra_keel_materials = None
            if rbm_mat is not None or shear_mat is not None or bearing_mat is not None or bolt_mat is not None:
                hydra_keel_materials = Data()
                if rbm_mat is not None:
                    hydra_keel_materials.uni_stress = rbm_mat.ultimate_tensile_strength
                    hydra_keel_materials.uni_rho = rbm_mat.density
                if shear_mat is not None:
                    hydra_keel_materials.bid_shear = shear_mat.ultimate_shear_strength
                    hydra_keel_materials.bid_rho = shear_mat.density
                if bearing_mat is not None:
                    hydra_keel_materials.bid_bearing = bearing_mat.ultimate_bearing_strength
                if bolt_mat is not None:
                    hydra_keel_materials.steel_shear = bolt_mat.ultimate_shear_strength

            fus = Hydra.compute_fuselage_weight(
                vehicle_mtow=MTOW, wing_span=maxSpan, min_wing_lift_fraction=min_wing_lift_fraction,
                fuselage_length=fuse.lengths.total, tech_factor=1.0,
                keel_materials=hydra_keel_materials,
            )
            fuse.mass_properties.center_of_gravity[0][0] = .45 * fuse.lengths.total
            fuse.mass_properties.mass = fus['total']
            weight.fuselage += fus['total']

        # -------------------------------------------------------------------------------
        # Boom (no Hydra model — Vahana/RCAIDE only)
        # -------------------------------------------------------------------------------
        for boom in vehicle.booms:
            boom_weight = Vahana.compute_boom_weight(boom) * Units.kg
            weight.booms += boom_weight
            boom.mass_properties.mass = boom_weight

        # -------------------------------------------------------------------------------
        # Pack up outputs — same schema as Physics_Based/NDARC
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
        # 2026-08-10: was hardcoded 0.0 - wired to n_pilots computed in the Payload section above.
        # getattr fallback: see note by flight_controls_constant above (settings here is
        # RCAIDE's own bare Data(), not WeightSettingsRCAIDEHydra).
        pilot_weight = getattr(settings, 'pilot_weight', 70.)
        output.operational_items.flight_crew = n_pilots * pilot_weight * Units.kg
        output.operational_items.flight_attendants = 0.0
        output.operational_items.passenger_service = 0.0
        output.operational_items.total = (output.operational_items.misc + output.operational_items.flight_crew
                                           + output.operational_items.flight_attendants + output.operational_items.passenger_service)

        # Contingency factor, 2026-08-10: applied to structure + systems + everything in
        # propulsion EXCEPT battery, matching AircraftDesign's own base/weight.py precedent
        # (AircraftWeight.calculate_weight(): `empty_weight *= contingency_factor` applied
        # BEFORE `energy_source` (battery) and passengers/pilot are added - user explicit
        # instruction "be consistent with aircraft design", not the initially-proposed
        # structure-only scope). operational_items.total (flight_crew/pilot) excluded too,
        # matching AircraftDesign's `weight.pilot` also being added after contingency.
        # getattr fallback: see note by flight_controls_constant above (settings here is
        # RCAIDE's own bare Data(), not WeightSettingsRCAIDEHydra).
        contingency_factor = getattr(settings, 'contingency_factor', 1.1)
        non_battery_empty = (output.empty.structural.total + output.empty.systems.total
                              + output.empty.propulsion.total - output.empty.propulsion.battery)
        non_battery_empty *= contingency_factor
        output.empty.total = non_battery_empty + output.empty.propulsion.battery + output.operational_items.total
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
