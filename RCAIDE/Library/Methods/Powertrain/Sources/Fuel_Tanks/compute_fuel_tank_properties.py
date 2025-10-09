# RCAIDE/Methods/Powertrain/Sources/Fuel_Tanks/compute_fuel_tank_properties.py
# 
# 
# Created:  Jul 2023, M. Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports
import RCAIDE
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  METHOD
# ----------------------------------------------------------------------------------------------------------------------  
def compute_fuel_tank_properties(tank,state,distributor,network_tag):
    '''
    UPDATE HEADER 
    ''' 
    
    if type(distributor) == RCAIDE.Library.Components.Powertrain.Distributors.Electrical_Bus:
        distributor_conditions = state.conditions.energy.busses[distributor.tag] 
    elif  type(distributor) == RCAIDE.Library.Components.Powertrain.Distributors.Fuel_Line: 
        distributor_conditions = state.conditions.energy.fuel_lines[distributor.tag]         
    
    tank_conditions = distributor_conditions.fuel_tanks[tank.tag]   
    mass_unknowns   = state.unknowns.network[tank.tag + '_mass'] 
     
    D   = state.numerics.time.differentiate      
    if len(D) > 0:
        tank_conditions.mass_flow_rate  =  distributor_conditions.fuel_mass_flow_rate * tank.flow_split_ratio
        dm_l                            = -tank_conditions.mass_flow_rate 
        R                               = np.dot(D,mass_unknowns)[:, 0] - dm_l[:, 0]
        R[0]                            = mass_unknowns[0] - tank_conditions.mass[0, 0]
    
        state.residuals.network[tank.tag + '_mass'] = R
        tank_conditions.mass[1:, 0]                 = mass_unknowns[1:, 0]

    return 