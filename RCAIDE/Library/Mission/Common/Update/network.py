# RCAIDE/Library/Missions/Common/Update/network.py
# 
# 
# Created:  Sep 2025, S Shekar

from re import S
from RCAIDE.Framework.Core import Data

from RCAIDE.Framework.Mission.Common import Conditions
import  numpy as np
from scipy.optimize import fsolve, least_squares

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
    
        unknown_keys = list(segment.state.unknowns.network[network.tag].keys()) 
        full_unkn_vals = Data()
        unknown_value  = Data()
        
        for unkn in unknown_keys:
            unknown_value[unkn]  = segment.state.unknowns.network[network.tag][unkn]  
            full_unkn_vals[unkn] = unknown_value[unkn] 
    
        if segment.state.network_numerics.solver.type  == 'least_squares':       
            result = least_squares(energy_model.evaluate, 
                        full_unkn_vals.pack_array(),
                        args=(segment,network),
                        method= segment.state.network_numerics.solver.method,
                        verbose = 2 if segment.state.network_numerics.solver.print_output is True else 0,
                        xtol=segment.state.network_numerics.solver.tolerance_solution,) 

            segment.state.network_numerics.solver.converged = result.success
            if result.success is False:
                print('The network solver fails with exit condition: ',result.status)
        elif segment.state.network_numerics.solver.type  == 'root_finder':
            result,_,ier,error_message = fsolve(energy_model.evaluate, 
                        full_unkn_vals.pack_array(),
                        args=(segment,network),
                        xtol=segment.state.network_numerics.solver.tolerance_solution,
                        maxfev = segment.state.numerics.solver.max_evaluations,
                        epsfcn = segment.state.numerics.solver.step_size,
                        full_output = 1) 

            segment.state.network_numerics.solver.converged = False if ier == 0 else True
            if ier == 0:
                print('The network solver fails with exit condition: ',error_message)
                
        else: # If a network solver is not needed it will unpack values from the missino solver
            energy_model.evaluate(full_unkn_vals.pack_array(), segment, network)