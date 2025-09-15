# RCAIDE/Library/Missions/Common/Update/thrust.py
# 
# 
# Created:  Jul 2023, M. Clarke
import  RCAIDE
from RCAIDE.Framework.Core import Data

import  numpy as np
from scipy.optimize import least_squares

# ----------------------------------------------------------------------------------------------------------------------
#  Update Thrust
# ---------------------------------------------------------------------------------------------------------------------- 
def network(segment):
    """ Updates the thrust vector of the vehicle from the propulsors 
        
        Assumptions:
        N/A
        
        Inputs:
            None 
                 
        Outputs: 
            None
      
        Properties Used:
        N/A
                    
    """ 

    # unpack
    energy_model = segment.analyses.energy
    


    unknown_keys = list(segment.state.unknowns.network.keys()) 
    full_unkn_vals = Data()
    unknown_value  = Data()
    
    for unkn in unknown_keys:
        unknown_value[unkn]  = segment.state.unknowns.network[unkn]  
        full_unkn_vals[unkn] = unknown_value[unkn] 

    initial_values    = full_unkn_vals.pack_array()        

    sol = least_squares(energy_model.evaluate, initial_values, args=([segment.state]),xtol=1e-14) 
    print(sol.x)
    a = 0

    # evaluate
    energy_model.evaluate(segment.state)