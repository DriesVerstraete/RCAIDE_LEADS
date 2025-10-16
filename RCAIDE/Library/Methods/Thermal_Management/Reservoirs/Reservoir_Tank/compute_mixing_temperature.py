# RCAIDE/Library/Methods/Thermal_Management/Reservoirs/Reservoir_Tank/compute_mixing_temperature.py


# Created:  Apr 2024, S. Shekar 

# ---------------------------------------------------------------------------------------------------------------------- 
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 
import  RCAIDE  

def compute_mixing_temperature(reservoir, state, coolant_line):
    """
    Computes the resultant temperature of the reservoir at each time step with coolant exchanging heat to the environment.

    :param reservoir: Reservoir Data Structure
        - reservoir.surface_area
        - reservoir.volume
        - reservoir.thickness
        - reservoir.material.conductivity
        - reservoir.material.emissivity
        - reservoir.coolant
    :type reservoir: dict
    :param state: State Data Structure
        - state.conditions.freestream.temperature
        - state.conditions.energy.coolant_line[reservoir.tag].coolant_temperature
    :type state: dict
    :param coolant_line: Coolant Line Data Structure
    :type coolant_line: dict 
    :type t_idx: int
    :return: Updated temperature of the reservoir coolant
    :rtype: float

    :Assumptions: 
        N/A

    :Source:
        None
    """  
    T_current = 0
    volume    = 0
    for reservoir in coolant_line.reservoirs:
        T_current += state.conditions.energy.coolant_lines[coolant_line.tag][reservoir.tag].coolant_temperature
        volume    += reservoir.volume
    
    T_current = T_current / len(coolant_line.reservoirs)

    # Reservoir Properties
    coolant       = reservoir.coolant
    rho_coolant   = coolant.compute_density(T_current)
    Cp_RES        = coolant.compute_cp(T_current)
    mass_coolant  = rho_coolant * volume 

    for battery in coolant_line.battery_modules: # THIS WILL BE REWRITTEN 
        for HAS in battery:
            if isinstance(HAS, RCAIDE.Library.Components.Thermal_Management.Batteries.Liquid_Cooled_Wavy_Channel):
                mass_flow_HAS = state.conditions.energy.coolant_lines[coolant_line.tag][HAS.tag].coolant_mass_flow_rate
                T_outlet_HAS  = state.conditions.energy.coolant_lines[coolant_line.tag][HAS.tag].outlet_coolant_temperature
                Cp_HAS        = coolant.compute_cp(T_outlet_HAS)

    for HEX in coolant_line.heat_exchangers:
        mass_flow_HEX = state.conditions.energy.coolant_lines[coolant_line.tag][HEX.tag].coolant_mass_flow_rate
        T_outlet_HEX  = state.conditions.energy.coolant_lines[coolant_line.tag][HEX.tag].outlet_coolant_temperature
        Cp_HEX        = coolant.compute_cp(T_outlet_HEX)

    # Solve for T_final using fsolve 
    # Ambient Air Temperature
    T_ambient = state.conditions.freestream.temperature

    # Compute heat loss to the environment
    # Properties of Reservoir
    A_surface       = reservoir.surface_area
    thickness       = reservoir.thickness
    conductivity    = reservoir.material.conductivity
    emissivity_res  = reservoir.material.emissivity

    # Heat Transfer properties
    sigma           = 5.69e-8  # Stefan Boltzmann Constant
    h               = 1000  # [W/m^2-K]
    emissivity_air  = 0.9

    # Heat Transfer due to conduction
    dQ_dt_cond = conductivity * A_surface * (T_final - T_ambient) / thickness

    # Heat Transfer due to natural convection
    dQ_dt_conv = h * A_surface * (T_final - T_ambient)

    # Heat Transfer due to radiation
    dQ_dt_rad = sigma * A_surface * ((emissivity_res * T_final ** 4) - (emissivity_air * T_ambient ** 4))

    dQ_dt_env = dQ_dt_cond + dQ_dt_conv + dQ_dt_rad  

    # Update the reservoir temperature
    state.conditions.energy.coolant_lines[coolant_line.tag][reservoir.tag].coolant_temperature  = T_final
    return 