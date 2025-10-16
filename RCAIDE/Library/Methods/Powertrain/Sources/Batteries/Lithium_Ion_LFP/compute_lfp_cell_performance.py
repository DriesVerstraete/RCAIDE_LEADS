# RCAIDE/Methods/Powertrain/Sources/Batteries/Lithium_Ion_LFP/compute_lfp_cell_performance.py
# 
# 
# Created: Nov 2024, S. Shekar

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
from RCAIDE.Framework.Core import Units
import numpy as np  
from copy import deepcopy

# ----------------------------------------------------------------------------------------------------------------------
# compute_lfp_cell_performance
# ----------------------------------------------------------------------------------------------------------------------  
def compute_lfp_cell_performance(battery_module, state, bus, coolant_lines):
    """
    Computes the performance of lithium iron phosphate (LFP) battery cells.
    
    Parameters
    ----------
    battery_module : BatteryModule
        The battery module containing LFP cells with the following attributes:
            - cell.electrode_area : float
                Area of the cell electrode [m²]
            - cell.surface_area : float
                Surface area of the cell [m²]
            - cell.mass : float
                Mass of the cell [kg]
            - cell.specific_heat_capacity : float
                Specific heat capacity of the cell [J/(kg·K)]
            - cell.discharge_performance_map : function
                Function that maps C-rate, temperature, and discharge capacity to voltage
            - electrical_configuration.series : int
                Number of cells in series
            - electrical_configuration.parallel : int
                Number of cells in parallel
    state : State
        Current system state containing conditions for all components
    bus : ElectricalBus
        The electrical bus connected to the battery module
    coolant_lines : list
        List of coolant lines that may be connected to the battery module 
    delta_t : numpy.ndarray
        Time step array [s]
    
    Returns
    -------
    stored_results_flag : bool
        Flag indicating if results were stored
    stored_battery_tag : str
        Tag of the battery module for which results were stored
    
    Notes
    -----
    This function computes the electrical and thermal performance of LFP battery cells
    within a battery module. It calculates:
        1. Current flow through the module and individual cells
        2. Heat generation due to joule heating and entropy changes
        3. Cell voltage under load based on state of charge, temperature, and current
        4. Power flow through the module
        5. Temperature changes based on heat generation and cooling (if present)
        6. State of charge and depth of discharge updates
        7. Charge throughput tracking
    
    **Major Assumptions**
        * All battery cells in a module exhibit the same thermal behavior
        * The cell temperature is assumed to be the temperature of the entire module
    
    **Theory**
    
    Heat generation in the cell includes both joule heating and entropic effects:
    
    .. math::
        \\dot{Q}_{joule} = \\frac{i_{cell}^2}{\sigma}
    
    .. math::
        \\dot{Q}_{entropy} = f(SOC) \\cdot T \\cdot \\frac{dU}{dT}
    
    where:
        - :math:`i_{cell}` is the current density [A/m²]
        - :math:`\\sigma` is the electrical conductivity [S/m]
        - :math:`SOC` is the state of charge
        - :math:`T` is the temperature [K]
        - :math:`\\frac{dU}{dT}` is the entropic coefficient
    
    See Also
    --------
    RCAIDE.Library.Components.Powertrain.Sources.Battery_Modules.Lithium_Ion_LFP
    """ 

    # ---------------------------------------------------------------------------------     
    # Time discretization
    # --------------------------------------------------------------------------------- 
    D = state.numerics.time.differentiate
    
    # ---------------------------------------------------------------------------------    
    # battery cell properties
    # --------------------------------------------------------------------------------- 
    electrode_area            = battery_module.cell.electrode_area 
    As_cell                   = battery_module.cell.surface_area
    cell_mass                 = battery_module.cell.mass    
    Cp                        = battery_module.cell.specific_heat_capacity       
    battery_module_data       = battery_module.cell.discharge_performance_map
    
    # ---------------------------------------------------------------------------------
    # Compute Bus electrical properties 
    # ---------------------------------------------------------------------------------    
    bus_conditions              = state.conditions.energy.busses[bus.tag]
    bus_config                  = bus.battery_module_electric_configuration 
    psi                         = state.conditions.energy.battery_fuel_cell_power_split_ratio
    P_bus                       = bus_conditions.power_draw*psi
    I_bus                       = bus_conditions.current_draw*psi   
    
    # ---------------------------------------------------------------------------------
    # Compute battery_module Conditions
    # -------------------------------------------------------------------------    
    battery_module_conditions = state.conditions.energy.busses[bus.tag].battery_modules[battery_module.tag]  
   
    E_module_max       = battery_module.maximum_energy * battery_module_conditions.cell.capacity_fade_factor 
    V_oc_module        = battery_module_conditions.voltage_open_circuit
    P_module           = battery_module_conditions.power
    P_cell             = battery_module_conditions.cell.power 
    Q_heat_module      = battery_module_conditions.heat_energy_generated
    Q_heat_cell        = battery_module_conditions.cell.heat_energy_generated 
    V_ul_cell          = battery_module_conditions.cell.voltage_under_load 
    I_module           = battery_module_conditions.current 
    I_cell             = battery_module_conditions.cell.current    

    # ---------------------------------------------------------------------------------                   
    # set unknowns 
    # ---------------------------------------------------------------------------------
    T_cell_unkn   = state.unknowns.network[battery_module.tag +  '_cell_temperature']
    SOC_cell_unkn = state.unknowns.network[battery_module.tag + '_cell_state_of_charge']
        
    
    # ---------------------------------------------------------------------------------
    # Compute battery_module electrical properties 
    # -------------------------------------------------------------------------    
    # Calculate the current going into one cell  
    n_series          = battery_module.electrical_configuration.series
    n_parallel        = battery_module.electrical_configuration.parallel 
    n_total           = n_series * n_parallel
    no_modules        = len(bus.battery_modules)

    # Scaling factors for numerical conditioning
    T_scale = 310.0
    E_scale = E_module_max
    
    # Scaled and bounded unknowns
    T_cell_scaled = T_cell_unkn / T_scale
    SOC_bounded   = np.clip(SOC_cell_unkn, 1e-4, 1.0)
        
    # ---------------------------------------------------------------------------------
    # Examine Thermal Management System
    # ---------------------------------------------------------------------------------
    HAS = None  
    for coolant_line in coolant_lines:
        for tag, item in  coolant_line.items():
            if tag == 'battery_modules':
                for sub_tag, sub_item in item.items():
                    if sub_tag == battery_module.tag:
                        for btms in  sub_item:
                            HAS = btms    


    # ---------------------------------------------------------------------------------------------------
    # Current State 
    # ---------------------------------------------------------------------------------------------------
    if bus_config == 'Series':
        I_module      = I_bus
    elif bus_config  == 'Parallel':
        I_module      = I_bus / len(bus.battery_modules)

       
    # ---------------------------------------------------------------------------------
    # Compute battery_module cell temperature 
    # ---------------------------------------------------------------------------------
    # Determine temperature increase         
    sigma           =  130  
    I_cell          = I_module / n_parallel   
    i_cell          = I_cell/electrode_area # current intensity (A/m²)
    q_dot_entropy   = (4.6810 * SOC_bounded**4 + (-8.3729) * SOC_bounded**3 + 3.7197 * SOC_bounded**2 + 0.4356 * SOC_bounded+ (-0.3027)) # Obtained from curve fitting the dUdt curve  
    q_dot_joule     = (i_cell**2)/(sigma)          
    Q_heat_cell     = (q_dot_joule + q_dot_entropy)*As_cell 
    Q_heat_module   = Q_heat_cell*n_total  
    V_ul_cell       = compute_lfp_cell_state(battery_module,battery_module_data,SOC_bounded,T_cell_unkn,abs(I_cell)) 
 
    # Power calculations
    P_module        = P_bus /no_modules + np.abs(Q_heat_module) 
    P_cell          = P_module/n_total
    
    if HAS is not None:
        dT_dt_scaled = HAS.compute_thermal_performance(battery_module, bus, coolant_line, Q_heat_cell,T_cell_unkn,T_scale,state)
    else:
        # Temperature residual with scaling
        dT_dt_scaled = Q_heat_cell / (cell_mass * Cp * T_scale)
    R_temp = np.dot(D, T_cell_scaled)[:, 0] - dT_dt_scaled[:, 0]
    R_temp[0] = T_cell_scaled[0] - battery_module_conditions.cell.temperature[0, 0] / T_scale
    state.residuals.network[battery_module.tag+ '_cell_temperature'] = R_temp
        
    # SOC residual with better conditioning
    dE_dt = -P_module
    R_soc = np.dot(D, SOC_cell_unkn * E_scale)[:, 0] - dE_dt[:, 0]
    R_soc[0] = SOC_cell_unkn[0] - battery_module_conditions.cell.state_of_charge[0, 0]
    state.residuals.network[battery_module.tag+ '_cell_state_of_charge'] = R_soc 
    
    # Update states
    battery_module_conditions.voltage_under_load            = V_ul_cell * n_series
    battery_module_conditions.cell.voltage_under_load       = V_ul_cell 
    battery_module_conditions.voltage_open_circuit          = V_oc_module 
    battery_module_conditions.cell.power                    = P_cell
    battery_module_conditions.cell.current                  = I_cell
    battery_module_conditions.current                       = I_module 
    battery_module_conditions.heat_energy_generated         = Q_heat_module
    battery_module_conditions.cell.heat_energy_generated    = Q_heat_cell 
    battery_module_conditions.cell.state_of_charge[1:,0]    = SOC_cell_unkn[1:,0]
    battery_module_conditions.state_of_charge[1:,0]         = SOC_cell_unkn[1:,0] 
    battery_module_conditions.cell.temperature[1:,0]        = T_cell_unkn[1:,0]
    battery_module_conditions.temperature[1:,0]             = T_cell_unkn[1:,0] 
    battery_module_conditions.cell.depth_of_discharge[1:,0] = 1. - SOC_cell_unkn[1:,0]
    battery_module_conditions.cell.energy[1:,0]             = SOC_cell_unkn[1:,0] * E_module_max / n_total
    battery_module_conditions.energy[1:,0]                  = SOC_cell_unkn[1:,0] * E_module_max
    
    # Charge throughput
    Q_prior = battery_module_conditions.cell.charge_throughput[0]
    dt      = np.diff(state.numerics.time.control_points[:,0])
    avg_I   = (I_cell[:-1, 0] + I_cell[1:, 0]) / 2
    Q_Ah    = np.atleast_2d(np.concatenate(([0.0], np.cumsum(dt*avg_I)))).T / Units.hr
    battery_module_conditions.cell.charge_throughput = Q_prior + Q_Ah
    
    stored_results_flag = True
    stored_battery_module_tag = battery_module.tag
    
    return stored_results_flag, stored_battery_module_tag 
 
 
def reuse_stored_lfp_cell_data(battery_module,state,bus,stored_results_flag, stored_battery_tag):
    """Reuses results from one propulsor for identical batteries       
    """
   
    state.conditions.energy.busses[bus.tag].battery_modules[battery_module.tag] = deepcopy(state.conditions.energy.busses[bus.tag].battery_modules[stored_battery_tag])      
    return


def compute_lfp_cell_state(battery_module, battery_module_data, SOC, T, I):
    """
    Computes the electrical state variables of a lithium iron phosphate (LFP) battery cell using look-up tables.
    
    Parameters
    ----------
    battery_module : BatteryModule
        The battery module containing LFP cells with the following attributes:
            - cell.nominal_capacity : float
                Nominal capacity of the cell [Ah]
    battery_module_data : function
        Look-up function that maps C-rate, temperature, and discharge capacity to voltage
    SOC : numpy.ndarray
        State of charge of the cell [unitless, 0-1]
    T : numpy.ndarray
        Battery cell temperature [K]
    I : numpy.ndarray
        Battery cell current [A]
    
    Returns
    -------
    V_ul : numpy.ndarray
        Under-load voltage [V]
    
    Notes
    -----
    This function computes the voltage of an LFP battery cell under load conditions
    by using a look-up table approach. It converts the state of charge to discharge
    capacity, calculates the C-rate, and then uses these values along with temperature
    to determine the cell voltage.
    
    The function applies limits to ensure the inputs are within the valid range of the
    look-up data:
        - SOC is limited to [0, 1]
        - Temperature is limited to [-10°C, 60°C]
        - Current is limited to [0A, 52A]
    
    **Major Assumptions**
        * The model is valid only within the specified temperature and current ranges
    
    See Also
    --------
    RCAIDE.Library.Components.Powertrain.Sources.Battery_Modules.Lithium_Ion_LFP
    """

    # Make sure things do not break by limiting current, temperature and current 
    capacity      = battery_module.cell.nominal_capacity
    SOC[SOC < 0.]   = 0.  
    SOC[SOC > 1.]   = 1.    
    DOD             = 1 - SOC 
    discharge_capacity = DOD*capacity
    

    T              = T-273
    # Operating Limits of the cell
    T[T<-10]       = -10 # model does not fit for below -10  degrees
    T[T>60]        =  60 # model does not fit for above 60 degrees

    
    I[I<0.0]      = 0.0
    I[I>52.0]     = 52.0
    C_rate        = I/capacity
     

    V_ul  = battery_module_data(C_rate, T, discharge_capacity)
    
    return V_ul