# RCAIDE/Methods/Powertrain/Sources/Fuel_Tanks/append_fuel_tank_conditions.py
# 
# 
# Created:  Jul 2023, M. Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports
import  RCAIDE
from RCAIDE.Framework.Core import Data
from RCAIDE.Framework.Mission.Common     import   Conditions, Residuals, Unknowns

# ----------------------------------------------------------------------------------------------------------------------
#  METHOD
# ----------------------------------------------------------------------------------------------------------------------  
def append_fuel_tank_conditions(tank, segment, distributor):
    """
    Appends initial conditions for fuel tank component during later mission analysis.
    
    Parameters
    ----------
    tank : FuelTank
        The fuel tank component for which conditions are being initialized.
    segment : Segment
        The mission segment in which the fuel tank is operating.
    distributor : Fuel Line of Bus
        The fuel line or bus connected to the fuel tank.
    
    Returns
    -------
    None
    
    Notes
    -----
    This function initializes the conditions for a fuel tank component at the start
    of a mission segment. It creates a Conditions object for the fuel tank within
    the segment's energy conditions dictionary, indexed by the fuel line tag and tank tag.
    
    The function initializes arrays for:
        - Mass flow rate [kg/s]
        - Mass [kg]
    
    These arrays are initialized with ones of the same length as the segment's state vector,
    which will be updated during mission analysis based on the fuel tank's performance
    and fuel consumption.
    
    See Also
    --------
    RCAIDE.Library.Components.Powertrain.Sources.Fuel_Tanks 
    """
    ones_row    = segment.state.ones_row
    
    if type(distributor) == RCAIDE.Library.Components.Powertrain.Distributors.Electrical_Bus: 
        distributor_conditions = segment.state.conditions.energy.busses[distributor.tag]
    elif  type(distributor) == RCAIDE.Library.Components.Powertrain.Distributors.Fuel_Line: 
        distributor_conditions = segment.state.conditions.energy.fuel_lines[distributor.tag]
        
    distributor_conditions.fuel_tanks[tank.tag]                           = Conditions()  
    distributor_conditions.fuel_tanks[tank.tag].mass                      = 0 * ones_row(1) 
    distributor_conditions.fuel_tanks[tank.tag].mass_flow_rate            = 0 * ones_row(1)  
    distributor_conditions.fuel_tanks[tank.tag].surface_temperature       = 0 * ones_row(1) 
    distributor_conditions.fuel_tanks[tank.tag].secondary_fuel_flow_rate  = tank.secondary_fuel_flow_rate * ones_row(1) 

         
    return 


def append_fuel_tank_segment_conditions(fuel_tank, segment, distributor): 

    if type(distributor) == RCAIDE.Library.Components.Powertrain.Distributors.Electrical_Bus: 
        distributor_conditions = segment.state.conditions.energy.busses[distributor.tag]
    elif  type(distributor) == RCAIDE.Library.Components.Powertrain.Distributors.Fuel_Line: 
        distributor_conditions = segment.state.conditions.energy.fuel_lines[distributor.tag]

    if segment.state.initials:  
        if type(distributor) == RCAIDE.Library.Components.Powertrain.Distributors.Electrical_Bus: 
            distributor_initals = segment.state.initals.conditions.energy.busses[distributor.tag]
        elif  type(distributor) == RCAIDE.Library.Components.Powertrain.Distributors.Fuel_Line: 
         distributor_initals = segment.state.initals.conditions.energy.fuel_lines[distributor.tag]
            
        distributor_conditions[fuel_tank.tag].mass[:,0]                     = distributor_initals[fuel_tank.tag].mass[-1,0]
    return

    
def append_fuel_tank_residual_and_unknowns(fuel_tank, segment, distributor,network):
    ones_row    = segment.state.ones_row

    distributor_unknowns   = segment.state.unknowns.network[network.tag].fuel_lines[distributor.tag]
    distributor_residuals = segment.state.residuals.network[network.tag].fuel_lines[distributor.tag]

    distributor_unknowns[fuel_tank.tag] = Unknowns()
    distributor_residuals[fuel_tank.tag] = Residuals()
    distributor_unknowns[fuel_tank.tag].mass  = ones_row(1) *fuel_tank.fuel.mass_properties.mass
    distributor_residuals[fuel_tank.tag].mass = ones_row(1)*0
    return
    
