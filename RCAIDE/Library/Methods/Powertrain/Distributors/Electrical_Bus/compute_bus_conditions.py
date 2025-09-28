# RCAIDE/Methods/Energy/Distributors/Electrical_Bus.py
# 
# 
# Created: Sep 2024, S. Shekar

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# imports 
import numpy as np

 
# ----------------------------------------------------------------------------------------------------------------------
# compute_bus_conditions
# ----------------------------------------------------------------------------------------------------------------------
def compute_bus_conditions(bus, state): 
    """
    Computes the conditions of the bus based on the response of the battery modules.
    
    Parameters
    ----------
    bus : ElectricalBus
        The electrical bus component with the following attributes:
            - tag : str
                Identifier for the bus
            - battery_modules : list
                List of battery modules connected to the bus
            - battery_module_electric_configuration : str
                Configuration of battery modules ('Series' or 'Parallel')
    state : State
        Current system state containing conditions for all components
    t_idx : int
        Current time index in the simulation
    delta_t : float
        Time step [s]
    
    Returns
    -------
    None
    
    Notes
    -----
    This function calculates the electrical properties of the bus based on the connected
    battery modules. It handles both series and parallel configurations of battery modules.
    
    For series configuration:
        - Voltage is the sum of all battery module voltages
        - Temperature and energy are averaged/summed across modules
        - State of charge is taken from the last battery module
    
    For parallel configuration:
        - Voltage is the same as the last battery module
        - Temperature and energy are averaged/summed across modules
        - State of charge is taken from the last battery module's cell
    
    The function also handles the fully charged state by setting charging current,
    power draw, and current draw to zero when the state of charge reaches 1.0.
    
    **Major Assumptions**
        * Battery modules are the primary energy source for the bus
        * In parallel configuration, all modules have the same voltage
    
    See Also
    --------
    RCAIDE.Library.Methods.Powertrain.Distributors.Electrical_Bus.append_bus_conditions
    RCAIDE.Library.Methods.Powertrain.Energy_Storage.Battery.compute_battery_module_conditions
    """
    bus_conditions = state.conditions.energy.busses[bus.tag]
    phi   = state.conditions.energy.hybrid_power_split_ratio
     
    if len(bus.battery_modules) != 0: 
        if bus.battery_module_electric_configuration == 'Series':
            bm_conditions                               = [bus_conditions.battery_modules[bm.tag] for bm in bus.battery_modules]
            bus_conditions.voltage_open_circuit   = sum(bm.voltage_open_circuit  for bm in bm_conditions)
            bus_conditions.voltage_under_load     = sum(bm.voltage_under_load  for bm in bm_conditions)
            bus_conditions.heat_energy_generated  = sum(bm.heat_energy_generated  for bm in bm_conditions)
            bus_conditions.efficiency             = (bus_conditions.power_draw *phi  + bus_conditions.heat_energy_generated )/(bus_conditions.power_draw *phi )
           
            bm_conditions                              = [bus_conditions.battery_modules[bm.tag] for bm in bus.battery_modules]
            bus_conditions.temperature         = sum(bm.temperature  for bm in bm_conditions)/ len(bus.battery_modules)
            bus_conditions.energy              = sum(bm.energy  for bm in bm_conditions)
            bus_conditions.state_of_charge     = bm_conditions[-1].state_of_charge 
    
        elif bus.battery_module_electric_configuration == 'Parallel':
            bm_conditions                               = [bus_conditions.battery_modules[bm.tag] for bm in bus.battery_modules]
            bus_conditions.heat_energy_generated  = sum(bm.heat_energy_generated  for bm in bm_conditions)
            bus_conditions.voltage_open_circuit   = bm_conditions[-1].voltage_open_circuit 
            bus_conditions.voltage_under_load     = bm_conditions[-1].voltage_under_load              
            bus_conditions.efficiency             = (bus_conditions.power_draw *phi  +  bus_conditions.heat_energy_generated )/(bus_conditions.power_draw *phi )

            bus_conditions.heat_energy_generated  = sum(bm.heat_energy_generated  for bm in bm_conditions)
            bus_conditions.temperature          = sum(bm.temperature  for bm in bm_conditions)/len(bus.battery_modules)
            bus_conditions.energy               = sum(bm.energy  for bm in bm_conditions)
            bus_conditions.state_of_charge      = bm_conditions[-1].cell.state_of_charge 
        

    # Handle fully charged state
    if state.conditions.energy.recharging and np.float16(bus_conditions.state_of_charge ) == 1:
        bus_conditions.charging_current  = 0
        bus_conditions.power_draw        = 0
        bus_conditions.current_draw      = 0
    return