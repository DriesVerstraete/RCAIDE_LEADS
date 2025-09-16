# RCAIDE/Library/Missions/Common/Update/network.py
# 
# 
# Created:  Sep 2025, S Shekar

from RCAIDE.Framework.Core import Data

import  numpy as np
from scipy.optimize import least_squares

# ----------------------------------------------------------------------------------------------------------------------
#  Solve Network
# ---------------------------------------------------------------------------------------------------------------------- 
def network(segment):
    """ Updates the **********
        
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
    for network  in segment.analyses.energy.vehicle.networks:
    
        unknown_keys = list(segment.state.unknowns[network.tag].keys()) 
        full_unkn_vals = Data()
        unknown_value  = Data()
        
        for unkn in unknown_keys:
            unknown_value[unkn]  = segment.state.unknowns[network.tag][unkn]  
            full_unkn_vals[unkn] = unknown_value[unkn] 

        if segment.state.network_numerics.solver.type  == 'least_squares':       
            result = least_squares(energy_model.evaluate, 
                        full_unkn_vals.pack_array(),
                        args=(segment.state,network),
                        method= segment.state.network_numerics.solver.method,
                        verbose = 2 if segment.state.network_numerics.solver.print_output is True else 0,
                        xtol=segment.state.network_numerics.solver.tolerance_solution,) 
            
            segment.state.network_numerics.solver.converged = result.success
            if result.success  is False:
                print(result.status)