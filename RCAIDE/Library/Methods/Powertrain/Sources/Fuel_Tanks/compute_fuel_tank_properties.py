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
    #####  distributor_conditions.fuel_tanks.fuel_tank.mass_flow_rate this is 0
    mass_unknowns =  state.unknowns[network_tag].fuel_lines[distributor.tag][tank.tag].mass[:,0]
    t0  = state.numerics.time.control_points[0][0]
    tf = state.numerics.time.control_points[-1][0]
    D = state.numerics.dimensionless.differentiate
    D_t = D/(tf - t0) 

    distributor_conditions.fuel_tanks.fuel_tank.mass_flow_rate  =  distributor_conditions.fuel_mass_flow_rate * tank.flow_split_ratio
    dm_l = -distributor_conditions.fuel_tanks.fuel_tank.mass_flow_rate 
    R = D_t @ mass_unknowns - dm_l[:,0]
    R[0] = mass_unknowns[0] - tank_conditions.mass[0][0]

    state.residuals[network_tag].fuel_lines[distributor.tag][tank.tag].mass = R
    tank_conditions.mass[1:,0] = mass_unknowns[1:]

    return 