# RCAIDE/Methods/Powertrain/Sources/Fuel_Tanks/append_fuel_tank_unknown_and_residual.py
# 
# 
# Created:  Nov 2025, M. Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  METHOD
# ----------------------------------------------------------------------------------------------------------------------  
def append_fuel_tank_unknown_and_residual(fuel_tank, segment):
    ones_row    = segment.state.ones_row
    segment.state.number_of_network_unknowns  += 1 
    segment.state.number_of_network_residuals += 1
    
    # unknown 
    segment.state.unknowns.network[fuel_tank.tag + '_mass']              = ones_row(1) *fuel_tank.fuel.mass_properties.mass 
    segment.state.residuals.network[fuel_tank.tag + '_mass']             = ones_row(1)*0 
    segment.state.unknowns_lower_bounds.network[fuel_tank.tag + '_mass'] = -np.inf * ones_row(1) 
    segment.state.unknowns_upper_bounds.network[fuel_tank.tag + '_mass'] = np.inf * ones_row(1)

    return
    
