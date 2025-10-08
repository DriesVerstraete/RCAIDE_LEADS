# RCAIDE/Methods/Powertrain/Sources/Fuel_Tanks/append_fuel_tank_unknown_and_residual.py
# 
# 
# Created:  Nov 2025, M. Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports
import  RCAIDE
from RCAIDE.Framework.Mission.Common     import   Conditions, Residuals, Unknowns

import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  METHOD
# ----------------------------------------------------------------------------------------------------------------------  
def append_fuel_tank_unknown_and_residual(fuel_tank, segment, distributor,network):
    ones_row    = segment.state.ones_row
    segment.state.number_of_network_unknowns  += 1 
    segment.state.number_of_network_residuals += 1
    
    # unknown 
    tank_unknowns   = segment.state.unknowns.network[network.tag].fuel_lines[distributor.tag].fuel_tanks  = Unknowns() 
    tank_residuals  = segment.state.residuals.network[network.tag].fuel_lines[distributor.tag].fuel_tanks = Residuals()
    distributor_lower_bounds   = segment.state.unknowns_lower_bounds.network[network.tag].fuel_lines[distributor.tag].fuel_tanks = Conditions()
    distributor_upper_bounds   = segment.state.unknowns_upper_bounds.network[network.tag].fuel_lines[distributor.tag].fuel_tanks = Conditions() 
    tank_unknowns[fuel_tank.tag].mass  = ones_row(1) *fuel_tank.fuel.mass_properties.mass
    
    # residual 
    tank_residuals[fuel_tank.tag].mass = ones_row(1)*0

    # lower bound 
    distributor_lower_bounds[fuel_tank.tag]      = Conditions()
    distributor_lower_bounds[fuel_tank.tag].mass = -np.inf * ones_row(1)
     
    # upper bound 
    distributor_upper_bounds[fuel_tank.tag]      = Conditions() 
    distributor_upper_bounds[fuel_tank.tag].mass = np.inf * ones_row(1)

    return
    
