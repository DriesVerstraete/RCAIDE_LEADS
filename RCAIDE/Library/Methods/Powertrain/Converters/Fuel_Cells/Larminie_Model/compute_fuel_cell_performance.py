# RCAIDE/Library/Methods/Powertrain/Converters/Fuel_Cells/Common/compute_fuel_cell_performance.py 
# 
# Created: Jan 2025, M. Clarke 

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------  
from RCAIDE.Framework.Core import Units 
import numpy as np 
import scipy as sp
# ----------------------------------------------------------------------
#  Larminie Model to Compute Fuel Cell Performance
# ---------------------------------------------------------------------- 
def compute_fuel_cell_performance(fuel_cell_stack, state, bus, coolant_lines):
    """
    Computes the performance of a fuel cell stack using the Larminie-Dicks model.
    
    Parameters
    ----------
    fuel_cell_stack : RCAIDE.Components.Energy.Converters.Fuel_Cell_Stack
        The fuel cell stack component containing cell properties and electrical configuration
    state : RCAIDE.Framework.Mission.Common.State
        Container for mission segment conditions
    bus : RCAIDE.Components.Energy.Distribution.Electric_Bus
        The electric bus to which the fuel cell stack is connected
    coolant_lines : list
        List of coolant line components for thermal management
    t_idx : int
        Current time index in the simulation
    delta_t : float
        Time step size [s]
         
    Returns
    -------
    stored_results_flag : bool
        Flag indicating that results have been stored for potential reuse
    stored_fuel_cell_stack_tag : str
        Tag identifier of the fuel cell stack with stored results
    
    Notes
    -----
    This function implements the Larminie-Dicks model to calculate fuel cell performance
    based on current operating conditions. It determines the optimal current density
    that matches the required power output, then calculates voltage, efficiency,
    and fuel consumption.
    
    The function handles both series and parallel electrical configurations for
    connecting the fuel cell stack to the electric bus.
    
    **Major Assumptions**
        * Uniform temperature distribution across all cells
        * No transient effects (steady-state operation at each time step)
        * Hydrogen is the only fuel considered
        * Ideal gas behavior
    
    **Theory**
    
    The Larminie-Dicks model calculates cell voltage as:
    
    .. math::
        V = E_0 - A\\ln(j) - Rj - m\\exp(nj)
    
    where:
        - E₀ is the open circuit voltage
        - A is the activation loss coefficient
        - R is the ohmic resistance
        - m and n are mass transport loss coefficients
        - j is the current density
    
    The efficiency is calculated as:
    
    .. math::
        \\eta = \\frac{V}{E_{ideal}}
    
    References
    ----------
    [1] Larminie, J., & Dicks, A. (2003). Fuel Cell Systems Explained (2nd ed.). John Wiley & Sons Ltd.
    
    See Also
    --------
    RCAIDE.Library.Methods.Powertrain.Converters.Fuel_Cells.Larminie_Model.compute_voltage
    RCAIDE.Library.Methods.Powertrain.Converters.Fuel_Cells.Larminie_Model.compute_power_difference
    """
    # ---------------------------------------------------------------------------------    
    # fuel cell stack properties 
    # --------------------------------------------------------------------------------- 
    fuel_cell         = fuel_cell_stack.fuel_cell 
    n_series          = fuel_cell_stack.electrical_configuration.series
    n_parallel        = fuel_cell_stack.electrical_configuration.parallel 
    bus_config        = bus.fuel_cell_stack_electric_configuration
    n_total           = n_series*n_parallel  
        
    # ---------------------------------------------------------------------------------
    # Compute Bus electrical properties   
    # ---------------------------------------------------------------------------------
    bus_conditions              = state.conditions.energy.busses[bus.tag]
    fuel_cell_stack_conditions  = bus_conditions.fuel_cell_stacks[fuel_cell_stack.tag]
    phi                         = state.conditions.energy.hybrid_power_split_ratio 
    P_bus                       = bus_conditions.power_draw*phi     
    P_stack                     = P_bus /len(bus.fuel_cell_stacks) 
    P_cell                      = P_stack/ n_total  

    # ---------------------------------------------------------------------------------
    # Compute fuel cell performance  
    # ---------------------------------------------------------------------------------  
    current_density_ukn         = state.unknowns.network[ fuel_cell_stack.tag + '_current_density'] 
    P_cell_cal                  = compute_power(current_density_ukn, fuel_cell)       
    state.residuals.network[ fuel_cell_stack.tag + '_power']  =  P_cell - P_cell_cal 
     
    V_fuel_cell                 = compute_voltage(fuel_cell,current_density_ukn)    
    efficiency                  = np.divide(V_fuel_cell, fuel_cell.ideal_voltage)
    mdot_cell                   = np.divide(P_cell,np.multiply(fuel_cell.propellant.specific_energy,efficiency)) 
    
    I_cell = P_cell / V_fuel_cell
    I_stack = I_cell * n_parallel
    if bus_config == 'Series':
        bus_conditions.current_draw = I_stack  
    elif bus_config  == 'Parallel': 
        bus_conditions.current_draw = I_stack * len(bus.fuel_cell_stacks)  
    
    fuel_cell_stack_conditions.power                                = P_stack
    fuel_cell_stack_conditions.current                              = I_stack
    fuel_cell_stack_conditions.voltage_open_circuit                 = V_fuel_cell *  n_series # assumes no losses
    fuel_cell_stack_conditions.voltage_under_load                   = V_fuel_cell *  n_series
    fuel_cell_stack_conditions.fuel_cell.voltage_open_circuit       = V_fuel_cell   # assumes no losses
    fuel_cell_stack_conditions.fuel_cell.voltage_under_load         = V_fuel_cell
    fuel_cell_stack_conditions.fuel_cell.power                      = P_cell
    fuel_cell_stack_conditions.fuel_cell.current                    = P_cell / V_fuel_cell 
    fuel_cell_stack_conditions.fuel_cell.inlet_H2_mass_flow_rate    = mdot_cell  
    fuel_cell_stack_conditions.H2_mass_flow_rate                    = mdot_cell * n_total # add fuel Line tag 
    
    stored_results_flag            = True
    stored_fuel_cell_stack_tag     = fuel_cell_stack.tag  

    return  stored_results_flag, stored_fuel_cell_stack_tag
 

def compute_power(current_density, fuel_cell_stack,sign=1.0):
    '''
    Function that determines the power output per cell, based on in 
    input current density
    
    Assumptions:
    None(calls other functions)
    
    Inputs:
    current_density      [Amps/m**2]
    fuel cell.
        interface area   [m**2]
        
    Outputs:
    power_out            [W]
    
    '''
    
    # sign variable is used so that you can maximize the power, by minimizing the -power
    i1            = current_density
    A             = fuel_cell_stack.fuel_cell.interface_area
    v             = compute_voltage(fuel_cell_stack,current_density)   # useful voltage vector
    power_out     = sign* np.multiply(v,i1)*A                    # obtain power output in W/cell 
    return power_out

def compute_voltage(fuel_cell_stack, current_density):
    """
    Calculates fuel cell voltage based on current density using the Larminie-Dicks model.
    
    Parameters
    ----------
    fuel_cell : RCAIDE.Components.Energy.Converters.Fuel_Cell
        The fuel cell component containing electrochemical parameters
            - r : float
                Area-specific resistance [Ohms*cm²]
            - A1 : float
                Tafel slope [V]
            - m : float
                Mass transport loss coefficient [V]
            - n : float
                Mass transport loss exponential coefficient [cm²/A]
            - Eoc : float
                Open circuit voltage [V]
    current_density : float or array
        Current density [A/m²]
        
    Returns
    -------
    v : float or array
        Cell voltage [V]
    
    Notes
    -----
    This function implements the Larminie-Dicks semi-empirical model to calculate
    fuel cell voltage as a function of current density. The model accounts for
    activation losses, ohmic losses, and concentration losses.
    
    **Major Assumptions**
        * Voltage curve follows the Larminie-Dicks model form
        * Steady-state operation (no transient effects)
        * Uniform current distribution across the cell
        * Constant temperature operation
    
    **Theory**
    
    The Larminie-Dicks model calculates cell voltage as:
    
    .. math::
        V = E_{oc} - r \\cdot i - A_1 \\ln(i) - m \\exp(n \\cdot i)
    
    where:
        - :math:`E_{oc}` is the open circuit voltage
        - r is the area-specific resistance
        - A_1 is the Tafel slope for activation losses
        - m and n are parameters for mass transport losses
        - i is the current density
    
    References
    ----------
    [1] Larminie, J., & Dicks, A. (2003). Fuel Cell Systems Explained (2nd ed.). John Wiley & Sons Ltd.
    
    See Also
    --------
    RCAIDE.Library.Methods.Powertrain.Converters.Fuel_Cells.Larminie_Model.compute_power
    """
    r   = fuel_cell_stack.fuel_cell.r/(1000*(Units.cm**2))
    Eoc = fuel_cell_stack.fuel_cell.Eoc 
    A1  = fuel_cell_stack.fuel_cell.A1  
    m   = fuel_cell_stack.fuel_cell.m   
    n   = fuel_cell_stack.fuel_cell.n   
    
    i1 = current_density/(0.001/(Units.cm**2.)) # current density(mA cm^-2)
    v  = Eoc-r*i1-A1*np.log(i1)-m*np.exp(n*i1)  #useful voltage vector

    return v

def compute_maximum_current_density(fuel_cell_stack):
    lb   = 0.0001/(Units.cm**2.)    #lower bound on fuel cell current density
    ub   = 1.2/(Units.cm**2.)
    sign = -1. # used to minimize -power 
    max_current_density = sp.optimize.fminbound(objective, lb, ub, args=(fuel_cell_stack, sign)) 
    return max_current_density


def objective(i, fuel_cell_stack,sign): 
    Power =  compute_power(i,fuel_cell_stack, sign)    
    return Power
