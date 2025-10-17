# RCAIDE/Methods/Powertrain/Converters/Fuel_Cells/Larminie_Model/append_fuel_cell_unknown_and_residual.py
# 
# 
# Created:  Oct 2025, M. Clarke 

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports
from RCAIDE.Framework.Core import  Units 

# ----------------------------------------------------------------------------------------------------------------------
#  append_fuel_cell_unknown_and_residual
# ----------------------------------------------------------------------------------------------------------------------
def append_fuel_cell_unknown_and_residual(fuel_cell,segment):

    ones_row    = segment.state.ones_row   

    # Conditions for recharging fuel_cell   
    segment.state.unknowns.network[ fuel_cell.tag + '_current_density']                =  1.0* ones_row(1)   
    segment.state.residuals.network[ fuel_cell.tag + '_power']                         =  0* ones_row(1)  
    segment.state.unknowns_upper_bounds.network[ fuel_cell.tag + '_current_density']   = 1.2/(Units.cm**2.)  * ones_row(1) 
    segment.state.unknowns_lower_bounds.network[ fuel_cell.tag + '_current_density']   = 0.0001/(Units.cm**2.)* ones_row(1)  
    segment.state.number_of_network_unknowns  += 1
    segment.state.number_of_network_residuals += 1          
      
    
    return 
