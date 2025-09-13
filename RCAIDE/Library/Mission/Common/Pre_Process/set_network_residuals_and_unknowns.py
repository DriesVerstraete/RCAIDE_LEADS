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
        networks = segment.analyses.energy.vehicle.netowkrks
        ones_row    = segment.state.ones_row 

        for network in networks:
            # Lets do only Fuel Tanks for now will figure out the rest later om 
            for fuel_line in network.fuel_lines:        
                for fuel_tank in fuel_line.fuel_tanks:
                    fuel_tank.append_residual_and_unkowns






        #  segment.state.unknowns.bank_angle = ones_row(1) * ctrls.bank_angle.initial_guess_values[0][0]