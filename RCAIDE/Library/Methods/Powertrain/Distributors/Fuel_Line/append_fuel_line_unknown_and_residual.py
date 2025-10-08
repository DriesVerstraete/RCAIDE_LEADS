#  RCAIDE/Methods/Energy/Distributors/Fuel_Line/append_fuel_line_unknown_and_residual.py
# 
# Created: Sep 2024, S. Shekar

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports   
from RCAIDE.Framework.Mission.Common     import   Conditions, Residuals, Unknowns

# ----------------------------------------------------------------------------------------------------------------------
#  METHODS
# ---------------------------------------------------------------------------------------------------------------------- 
def append_fuel_line_unknown_and_residual(fuel_line,segment,network): 
    """ 
    """ 
    # append fuel line unknown and residual
    
    # for fuel tanks appended on fuel line     
    segment.state.unknowns.network[network.tag].fuel_lines[fuel_line.tag].fuel_tanks              = Unknowns() 
    segment.state.residuals.network[network.tag].fuel_lines[fuel_line.tag].fuel_tanks             = Residuals()
    segment.state.unknowns_lower_bounds.network[network.tag].fuel_lines[fuel_line.tag].fuel_tanks = Conditions()
    segment.state.unknowns_upper_bounds.network[network.tag].fuel_lines[fuel_line.tag].fuel_tanks = Conditions()     

    return