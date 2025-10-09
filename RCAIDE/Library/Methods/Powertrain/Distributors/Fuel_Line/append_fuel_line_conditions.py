#  RCAIDE/Methods/Energy/Distributors/Fuel_Line/append_fuel_line_conditions.py
# 
# Created: Sep 2024, S. Shekar

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports   
from RCAIDE.Framework.Mission.Common     import   Conditions
# ----------------------------------------------------------------------------------------------------------------------
#  METHODS
# ---------------------------------------------------------------------------------------------------------------------- 
def append_fuel_line_conditions(fuel_line,segment,network): 
    """
    Appends conditions for the fuel line to the segment's energy conditions dictionary.

    Parameters
    ----------
    fuel_line : RCAIDE.Library.Components.Distributors.Fuel_Line
    
    Returns
    -------
    None
        This function modifies the segment.state.conditions.energy dictionary in-place.
    
    Notes
    -----
    This function creates a Conditions object for the fuel line within the segment's
    energy conditions dictionary, indexed by the fuel line tag. It initializes various fuel line
    properties as zero arrays with the same length as the segment's state vector.
    
    The initialized properties include: 
        - Heat energy generated
        - Efficiency
        - Temperature
        - Energy
        - Flow Rate 
        - Regenerative power
    
    For segments with an initial battery state of charge specified, the function also
    sets the initial energy and state of charge values accordingly.
    
    See Also
    --------
    RCAIDE.Library.Methods.Powertrain.Distributors.Electrical_fuel_line.compute_fuel_line_conditions
    """
    ones_row                                                                     = segment.state.ones_row

    # ------------------------------------------------------------------------------------------------------            
    # Create fuel_line results data structure  
    # ------------------------------------------------------------------------------------------------------ 
    segment.state.conditions.energy.fuel_lines[fuel_line.tag]                                     = Conditions() 
    segment.state.conditions.energy.fuel_lines[fuel_line.tag].power_draw                          = 0 * ones_row(1)
    segment.state.conditions.energy.fuel_lines[fuel_line.tag].hybrid_power_split_ratio            = segment.hybrid_power_split_ratio * ones_row(1)  
    segment.state.conditions.energy.fuel_lines[fuel_line.tag].heat_energy_generated               = 0 * ones_row(1) 
    segment.state.conditions.energy.fuel_lines[fuel_line.tag].efficiency                          = 0 * ones_row(1)
    segment.state.conditions.energy.fuel_lines[fuel_line.tag].temperature                         = 0 * ones_row(1)
    segment.state.conditions.energy.fuel_lines[fuel_line.tag].energy                              = 0 * ones_row(1)  
    segment.state.conditions.energy.fuel_lines[fuel_line.tag].fuel_mass_flow_rate                 = 0 * ones_row(1)  
    segment.state.conditions.energy.fuel_lines[fuel_line.tag].fuel_tanks                          = Conditions() 


    
    m_total = 0
    for fuel_tank in fuel_line.fuel_tanks:
        m_total += fuel_tank.fuel.mass_properties.mass
    
    for fuel_tank in fuel_line.fuel_tanks: 
        if fuel_tank.fuel_selector_valve.fuel_flow_split_ratio == 0: 
            fuel_tank.fuel_selector_valve.fuel_flow_split_ratio = (fuel_tank.fuel.mass_properties.mass / m_total)  
 
    return