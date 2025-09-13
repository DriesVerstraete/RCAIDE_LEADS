# RCAIDE/Library/Missions/Common/Pre_Process/set_residuals_and_unknowns.py
# 
# 
# Created:  Jul 2023, M. Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------  
from RCAIDE.Framework.Core import Units

# ----------------------------------------------------------------------------------------------------------------------
#  set_residuals_and_unknowns
# ----------------------------------------------------------------------------------------------------------------------  
def set_network_residuals_and_unknowns(mission):
     
     
    
    for segment in mission.segments: 
        segment.state.number_of_network_unknowns = 0 
        networks = segment.analyses.energy.vehicle.networks
        for network in networks:
            # Lets do only Fuel Tanks for now will figure out the rest later om 
            for fuel_line in network.fuel_lines:        
                for fuel_tank in fuel_line.fuel_tanks:
                    segment.state.number_of_network_unknowns  += 1 
                    fuel_tank.append_residual_and_unkowns