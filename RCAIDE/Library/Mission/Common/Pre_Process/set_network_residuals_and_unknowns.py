# RCAIDE/Library/Missions/Common/Pre_Process/set_residuals_and_unknowns.py
# 
# 
# Created:  Jul 2023, M. Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  set_residuals_and_unknowns
# ----------------------------------------------------------------------------------------------------------------------  
def set_network_residuals_and_unknowns(mission):
     
    for segment in mission.segments: 
        segment.state.number_of_network_unknowns   = 0 
        segment.state.number_of_network_residuals  = 0
        networks = segment.analyses.energy.vehicle.networks
        for network in networks:                        

            # ---------------------------------------------------------------------------------------------
            # Propulsors 
            # ---------------------------------------------------------------------------------------------
            for p_i,propulsor in  enumerate(network.propulsors):  
                if network.identical_propulsors == True and p_i ==0:
                    propulsor.append_unknowns_and_residuals(segment)
                    
            # ---------------------------------------------------------------------------------------------            
            # Distributors 
            # ---------------------------------------------------------------------------------------------
            # Fuel Line 
            for fuel_line in network.fuel_lines:                 
                fuel_line.append_unknowns_and_residuals(segment,network)                
         
                for fuel_tank in fuel_line.fuel_tanks:
                    fuel_tank.append_unknowns_and_residuals(segment)                        
                 
            # Bus 
            for bus in network.busses:    
                bus.append_unknowns_and_residuals(segment,network)                  

                for bat_i,battery_module in  enumerate(bus.battery_modules):  
                    if bus.identical_battery_modules == True and bat_i ==0:
                        battery_module.append_unknowns_and_residuals(segment) 
        
                for fc_i,fuel_cell_stack in  enumerate(bus.fuel_cell_stacks):    
                    if bus.identical_fuel_cell_stacks == True and fc_i ==0:
                        fuel_cell_stack.append_unknowns_and_residuals(segment)    
    
            # # Ensure the mission knows how to pack and unpack the unknowns and residuals
            segment.process.iterate.unknowns.mission.network   = network.unpack_unknowns 
            segment.process.iterate.residuals.mission.network  = network.residuals