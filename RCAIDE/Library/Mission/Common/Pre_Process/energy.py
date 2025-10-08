# RCAIDE/Library/Missions/Common/Pre_Process/energy.py
# 
# 
# Created:  Jul 2023, M. Clarke
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE
# ----------------------------------------------------------------------------------------------------------------------  
import RCAIDE 
# ----------------------------------------------------------------------------------------------------------------------
#  energy
# ----------------------------------------------------------------------------------------------------------------------  
def energy(mission):
    """ Pre-processes energy network by appending all unknowns and residuals             
    """       
    for segment in mission.segments: 
        for network in segment.analyses.energy.vehicle.networks: 
            if type(network) == RCAIDE.Framework.Networks.Hybrid:
                if segment.hybrid_power_split_ratio == None:
                    raise AssertionError('Hybridization power split ratio not set! Specify in mission segment') 
                if segment.battery_fuel_cell_power_split_ratio == None:
                    raise AssertionError('Battery/Fuel cell power split ratio not set! Specify in mission segment')                 
            elif type(network) == RCAIDE.Framework.Networks.Fuel: 
                if segment.hybrid_power_split_ratio == None:                
                    segment.hybrid_power_split_ratio = 0.0
                    segment.battery_fuel_cell_power_split_ratio = 0.0
            elif type(network) == RCAIDE.Framework.Networks.Electric: 
                if segment.hybrid_power_split_ratio == None:                
                    segment.hybrid_power_split_ratio = 1.0  
                    segment.battery_fuel_cell_power_split_ratio = 1.0
            elif type(network) == RCAIDE.Framework.Networks.Fuel_Cell: 
                if segment.hybrid_power_split_ratio == None:                
                    segment.hybrid_power_split_ratio = 1.0  
                    segment.battery_fuel_cell_power_split_ratio = 0.0
            segment.state.conditions.energy.hybrid_power_split_ratio            = segment.hybrid_power_split_ratio * segment.state.ones_row(1)  
            segment.state.conditions.energy.battery_fuel_cell_power_split_ratio = segment.battery_fuel_cell_power_split_ratio * segment.state.ones_row(1)                    
            
        
            # ---------------------------------------------------------------------------------------------
            # Propulsors 
            # ---------------------------------------------------------------------------------------------
            for p_i,propulsor in  enumerate(network.propulsors): 
                propulsor.append_operating_conditions(segment)
        
            # ---------------------------------------------------------------------------------------------
            # Converters 
            # ---------------------------------------------------------------------------------------------    
            for converter in network.converters: 
                converter.append_operating_conditions(segment,network)                           
        
            # ---------------------------------------------------------------------------------------------            
            # Distributors 
            # ---------------------------------------------------------------------------------------------
            # Fuel Line 
            for fuel_line in network.fuel_lines:          
                fuel_line.append_operating_conditions(segment,network) 
        
                for fuel_tank in fuel_line.fuel_tanks:
                    fuel_tank.append_operating_conditions(segment,fuel_line,network)             
        
            # Bus 
            for bus in network.busses:   
                bus.append_operating_conditions(segment,network) 
        
                for bat_i,battery_module in  enumerate(bus.battery_modules): 
                    battery_module.append_operating_conditions(segment,bus)  
        
                for fc_i,fuel_cell_stack in  enumerate(bus.fuel_cell_stacks): 
                    fuel_cell_stack.append_operating_conditions(segment,bus,network)                                
        
                for tag, bus_item in bus.items():  
                    if issubclass(type(bus_item), RCAIDE.Library.Components.Component):
                        bus_item.append_operating_conditions(segment,bus,network)
                        
            # Coolant Line 
            for coolant_line in network.coolant_lines:
                coolant_line.append_operating_conditions(segment,network)       
        
                for battery_module in coolant_line.battery_modules: 
                    for btms in battery_module:
                        btms.append_operating_conditions(segment,coolant_line,network)
        
                for heat_exchanger in coolant_line.heat_exchangers: 
                    heat_exchanger.append_operating_conditions(segment,coolant_line,network)
        
                for reservoir in coolant_line.reservoirs: 
                    reservoir.append_operating_conditions(segment,coolant_line,network)
                     