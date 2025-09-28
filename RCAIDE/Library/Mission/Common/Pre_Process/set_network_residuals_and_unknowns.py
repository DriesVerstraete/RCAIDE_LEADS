# RCAIDE/Library/Missions/Common/Pre_Process/set_residuals_and_unknowns.py
# 
# 
# Created:  Jul 2023, M. Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------  
import RCAIDE 
# ----------------------------------------------------------------------------------------------------------------------
#  set_residuals_and_unknowns
# ----------------------------------------------------------------------------------------------------------------------  
def set_network_residuals_and_unknowns(mission): 
    for segment in mission.segments: 
        segment.state.number_of_network_unknowns   = 0 
        segment.state.number_of_network_residuals  = 0
        networks = segment.analyses.energy.vehicle.networks
        for network in networks: 
            for fuel_line in network.fuel_lines:                  
                fuel_line.append_operating_conditions(segment,network)  
                for fuel_tank in fuel_line.fuel_tanks:
                    fuel_tank.append_operating_conditions(segment,fuel_line,network)
            for bus in network.busses:                    
                bus.append_operating_conditions(segment,network)  
                for index,battery_module in enumerate(bus.battery_modules): 
                    battery_module.append_operating_conditions(segment,bus) 
                    if bus.identical_battery_modules == True and index ==0:
                        battery_module.append_unknowns_residuals(segment,bus,network) 
                for tag, bus_item in bus.items():  
                    if issubclass(type(bus_item), RCAIDE.Library.Components.Component):
                        bus_item.append_operating_conditions(segment,bus)
    