# RCAIDE/Library/Missions/Common/Pre_Process/set_residuals_and_unknowns.py
# 
# 
# Created:  Jul 2023, M. Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------  
from RCAIDE.Framework.Core import Units, Data
from RCAIDE.Framework.Mission.Common import Conditions
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
                segment.state.unknowns.network                 = Conditions()
                segment.state.unknowns.network.fuel_lines      = Conditions()
                segment.state.residuals.network                = Conditions()
                segment.state.residuals.network.fuel_lines     = Conditions()
       
                for fuel_tank in fuel_line.fuel_tanks:
                    segment.state.number_of_network_unknowns  += 1 
                    fuel_tank.append_residual_and_unknowns(segment,fuel_line)
                    a=0