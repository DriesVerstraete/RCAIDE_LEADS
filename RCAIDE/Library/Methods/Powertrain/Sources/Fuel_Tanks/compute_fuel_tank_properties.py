# RCAIDE/Methods/Powertrain/Sources/Fuel_Tanks/compute_fuel_tank_properties.py
# 
# 
# Created:  Jul 2023, M. Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports
import RCAIDE
from RCAIDE.Library.Mission.Common.Unpack_Unknowns.energy import unknowns
import numpy as np
# ----------------------------------------------------------------------------------------------------------------------
#  METHOD
# ----------------------------------------------------------------------------------------------------------------------  
def compute_fuel_tank_properties(tank,state,distributor,network_tag):
    # '''
    # SAI HEADER
    # ''' 
    
    if type(distributor) == RCAIDE.Library.Components.Powertrain.Distributors.Electrical_Bus:
        distributor_conditions = state.conditions.energy.busses[distributor.tag] 
    elif  type(distributor) == RCAIDE.Library.Components.Powertrain.Distributors.Fuel_Line: 
        distributor_conditions = state.conditions.energy.fuel_lines[distributor.tag]         
    
    tank_conditions = distributor_conditions.fuel_tanks[tank.tag]      
   
    mass_unknowns =  state.unknowns.network[network_tag].fuel_lines[distributor.tag][tank.tag].mass[:,0]

    t0  = state.numerics.time.control_points[0][0]
    tf = state.numerics.time.control_points[-1][0]
    D = state.numerics.dimensionless.differentiate
    D_t = D/(tf - t0) 

    tank_conditions.mass_flow_rate  =  distributor_conditions.fuel_mass_flow_rate * tank.flow_split_ratio
    dm_l = -tank_conditions.mass_flow_rate[:, 0]
    R = D_t @ mass_unknowns - dm_l
    R[0] = mass_unknowns[0] - tank_conditions.mass[0, 0]

    state.residuals.network[network_tag].fuel_lines[distributor.tag][tank.tag].mass = R
    tank_conditions.mass[1:, 0] = mass_unknowns[1:]

    return 