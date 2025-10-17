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
def append_coolant_line_conditions(coolant_line,segment): 
    """
     
    """
    ones_row                                                                     = segment.state.ones_row

    # ------------------------------------------------------------------------------------------------------            
    # Create fuel_line results data structure  
    # ------------------------------------------------------------------------------------------------------ 
    segment.state.conditions.energy.coolant_lines[coolant_line.tag]                                     = Conditions()  
     

    return
