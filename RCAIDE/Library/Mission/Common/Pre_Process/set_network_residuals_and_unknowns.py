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
        segment.state.number_of_network_unknowns   = 0 
        segment.state.number_of_network_residuals  = 0
        networks = segment.analyses.energy.vehicle.networks
        for network in networks:
            segment.state.unknowns.network[network.tag]                 = Conditions()
            segment.state.residuals.network[network.tag]                = Conditions()
            segment.state.lower_bounds.network[network.tag]                 = Conditions()
            segment.state.upper_bounds.network[network.tag]                = Conditions()
            # Lets do only Fuel Tanks for now will figure out the rest later om 
            for fuel_line in network.fuel_lines:   
                segment.state.unknowns.network[network.tag].fuel_lines      = Conditions()
                segment.state.residuals.network[network.tag].fuel_lines     = Conditions()     
                segment.state.lower_bounds.network[network.tag].fuel_lines      = Conditions()
                segment.state.upper_bounds.network[network.tag].fuel_lines     = Conditions()               
                fuel_line.append_operating_conditions(segment,network) 

                for fuel_tank in fuel_line.fuel_tanks:
                    fuel_tank.append_operating_conditions(segment,fuel_line,network)