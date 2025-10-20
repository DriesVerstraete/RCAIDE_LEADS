# RCAIDE/Framework/Networks/Network.py 
#
# Created:  Mar 2025, M.Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  Imports
# ---------------------------------------------------------------------------------------------------------------------
# RCAIDE Imports
import  RCAIDE
from RCAIDE.Library.Mission.Common.Unpack_Unknowns.energy import unknowns
from RCAIDE.Library.Methods.Powertrain.Systems.compute_avionics_power_draw                import compute_avionics_power_draw
from RCAIDE.Library.Methods.Powertrain.Systems.compute_systems_power_draw                 import compute_systems_power_draw
from RCAIDE.Library.Methods.Powertrain.Converters.Motor.compute_motor_performance         import *
from RCAIDE.Library.Methods.Powertrain.Converters.Generator.compute_generator_performance import * 
from RCAIDE.Library.Components import Component

# python imports 
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  Network
# ---------------------------------------------------------------------------------------------------------------------- 
class Network(Component):  
    """ Generalized Hybrid Energy Network (powertrain) Class capable of creating all derivatives of hybrid
    networks, the conventional fuel network and the all-electric network.
    
                                            GENERIC NETWORK
          .........................................:..........................................                          
          :                        :                                :                        :
    .-------------.         .-------------.                 .-------------.           .-------------.                     
    | propulsor 1 |         | propulsor 2 |                 | propulsor 2 |           | propulsor 3 | 
    '-------------'         '-------------'                 '-------------'           '-------------'            
          ||                       ||                              ||                        ||                                  
          ||   .-------------.     ||                              ||  .-------------.       ||
          ||== | converter 1 |====== electric bus / fuel line =========| converter 2 |=======|| 
               '-------------'                                         '-------------'  
                           
    Attributes
    ----------
    tag : str
        Identifier for the network   
    
    Notes
    -----
    The evaluate function is broken into three sections: Section 1 computes all the forces and moments
    from propulsors regardless of if they are powered by fuel or an electrochemical energy storage system;
    Section 2 computees the perfomrance of any converters on the distrution lines, for example,
    turboshafts, motors, pumps etc; and Section 3 computes the thermal mangement of the system as
    well as energy consumtion of the powertrain. The state of storage devices such as covnentional fuel tanks,
    batteries are also updates. Propulsor groups can be "active" or "inactive" to simulate
    engine out conditions. Energy consumtion from avionics is also modeled 
    
    **Definitions** 
    'Propulsor Group'
        Any single or group of Components that work together to provide thrust.
    
    See Also
    --------
    RCAIDE.Library.Framework.Networks.Fuel
        Fuel network class 
    RCAIDE.Library.Framework.Networks.Fuel_Cell
        Fuel_Cell network class 
    RCAIDE.Library.Framework.Networks.Electric
        All-Electric network class  
    """      
    
    def __defaults__(self):
        """ This sets the default values for the network to function.
        """        
        self.tag                          = 'network'
        self.propulsors                   = Container()  
        self.busses                       = Container()
        self.coolant_lines                = Container()
        self.fuel_lines                   = Container()
        self.converters                   = Container()
        self.identical_propulsors         = True 
        self.reverse_thrust               = False
        self.wing_mounted                 = True   
        self.system_voltage               = None  
        
    # linking the different network components
    def evaluate(network,state,center_of_gravity):

        """ Computes the performance of the network
        """  
        # unpack   
        conditions           = state.conditions 
        busses               = network.busses 
        fuel_lines           = network.fuel_lines 
        coolant_lines        = network.coolant_lines
        converters           = network.converters 
        total_thrust         = 0. * state.ones_row(3) 
        total_mech_power     = 0. * state.ones_row(1) 
        total_elec_power     = 0. * state.ones_row(1) 
        total_moment         = 0. * state.ones_row(3)  
        total_mdot           = 0. * state.ones_row(1)   
        reverse_thrust       = network.reverse_thrust 
    
        # ----------------------------------------------------------       
        # Section 1.0 Propulsor Performance 
        # ----------------------------------------------------------
        # 1.1 Fuel Propulsors  
        for fuel_line in fuel_lines: 
            conditions.energy.fuel_lines[fuel_line.tag].fuel_mass_flow_rate = 00
            for propulsor_group in fuel_line.assigned_propulsors:
                stored_results_flag  = False
                stored_propulsor_tag = None 
                conditions.energy.fuel_lines[fuel_line.tag].fuel_flow_rate = 0
                for propulsor_tag in propulsor_group:
                    propulsor = network.propulsors[propulsor_tag]
                    if propulsor.active and fuel_line.active:   
                        if network.identical_propulsors == False or stored_results_flag == False:
                            # run analysis  
                            T,M,P,P_elec,stored_results_flag,stored_propulsor_tag = propulsor.compute_performance(state, center_of_gravity= center_of_gravity)
                        else:              
                            # use previous propulsor results 
                            T,M,P,P_elec = propulsor.reuse_stored_data(state,network,stored_propulsor_tag=stored_propulsor_tag,center_of_gravity= center_of_gravity)
        
                        total_thrust      += T   
                        total_moment      += M   
                        total_mech_power  += P   
        
                        # compute total mass flow rate
                        conditions.energy.fuel_lines[fuel_line.tag].fuel_mass_flow_rate += conditions.energy.propulsors[propulsor.tag].fuel_mass_flow_rate
                
        # 1.2 Electric Propulsors         
        for bus in busses:            
            avionics             = bus.avionics 
            systems              = bus.systems 
    
            # Avionics Power Consumtion 
            compute_avionics_power_draw(avionics,bus,conditions) 
            compute_systems_power_draw(systems,bus,conditions) 
    
            # Bus Voltage 
            bus_voltage = bus.voltage * state.ones_row(1)       
    
            if conditions.energy.recharging:             
                bus.charging_current         = bus.nominal_capacity * bus.charging_c_rate 
                charging_power               = (bus.charging_current*bus_voltage*bus.power_split_ratio) 
                conditions.energy.busses[bus.tag].power_draw   -= charging_power/bus.efficiency
                conditions.energy.busses[bus.tag].current_draw  = -conditions.energy.busses[bus.tag].power_draw/bus.voltage
    
            else:
                for propulsor_group in bus.assigned_propulsors:
                    stored_results_flag  = False  
                    for propulsor_tag in propulsor_group:
                        propulsor =  network.propulsors[propulsor_tag]
                        if propulsor.active and bus.active:      
                            if network.identical_propulsors == False or stored_results_flag == False: 
                                # run analysis  
                                T,M,P_mech,P_elec,stored_results_flag,stored_propulsor_tag = propulsor.compute_performance(state,center_of_gravity= center_of_gravity)
                            else:    
                                T,M,P_mech,P_elec  = propulsor.reuse_stored_data(state,network,stored_propulsor_tag=stored_propulsor_tag,center_of_gravity=center_of_gravity)
    
                            total_thrust      += T   
                            total_moment      += M   
                            total_mech_power  += P_mech 
                            total_elec_power  += P_elec 
    
                # compute power from each componemnt 
                conditions.energy.busses[bus.tag].power_draw        += (total_elec_power- state.conditions.energy.busses[bus.tag].regenerative_power*bus_voltage ) * bus.power_split_ratio  /bus.efficiency   
                conditions.energy.busses[bus.tag].current_draw       = conditions.energy.busses[bus.tag].power_draw/bus_voltage  
             
        # ------------------------------------------------------------------------------------------------------------------- 
        # Section 2.0 Converters
        # -------------------------------------------------------------------------------------------------------------------  
        # 2.1 Fuel Converters         
        for fuel_line in fuel_lines: 
            for converter_group in fuel_line.assigned_converters:
                stored_conveter_tag = False
                for converter_tag in converter_group:
                    converter =  converters[converter_tag]
                    if converter.active and fuel_line.active: 
                        converter.inverse_calculation = True 
                        if isinstance(converter,RCAIDE.Library.Components.Powertrain.Converters.Turboelectric_Generator): 
                            if stored_conveter_tag is False:
                                generator             = converter.generator   
                                state.conditions.energy.converters[generator.tag].outputs.power  =  total_elec_power*(1 - state.conditions.energy.hybrid_power_split_ratio ) 
                                P_mech, P_elec, stored_results_flag,stored_conveter_tag          = converter.compute_performance(state,fuel_line,bus)  
                                conditions.energy.busses[bus.tag].power_draw                    -= P_elec/bus.efficiency
                                conditions.energy.fuel_lines[fuel_line.tag].fuel_mass_flow_rate += conditions.energy.converters[converter.tag].fuel_mass_flow_rate   
                            else:
                                generator             = converter.generator   
                                state.conditions.energy.converters[generator.tag].outputs.power  =  total_elec_power*(1 - state.conditions.energy.hybrid_power_split_ratio ) 
                                P_mech, P_elec                                                   = converter.reuse_stored_data(state,network,stored_conveter_tag,fuel_line,bus)  
                                conditions.energy.busses[bus.tag].power_draw                     -= P_elec/bus.efficiency
                                conditions.energy.fuel_lines[fuel_line.tag].fuel_mass_flow_rate  += conditions.energy.converters[converter.tag].fuel_mass_flow_rate   

                        if isinstance(converter,RCAIDE.Library.Components.Powertrain.Converters.Turboshaft):   
                            state.conditions.energy.converters[converter.tag].power     = total_mech_power*(1 - state.conditions.energy.hybrid_power_split_ratio )   
                            P_mech, P_elec,stored_results_flag,stored_propulsor_tag     = converter.compute_performance(state)   
                            conditions.energy.fuel_lines[fuel_line.tag].fuel_mass_flow_rate  += conditions.energy.converters[converter.tag].fuel_mass_flow_rate  
                    
        # 2.1 Electric Converters                            
        for bus in busses:         
            for converter_group in bus.assigned_converters:
                for converter_tag in converter_group:
                    converter =  converters[converter_tag]
                    if converter.active: 
                        converter.inverse_calculation = True
                        if isinstance(converter,RCAIDE.Library.Components.Powertrain.Converters.DC_Motor) or isinstance(converter,RCAIDE.Library.Components.Powertrain.Converters.PMSM_Motor):  
                            compute_motor_performance(converter,conditions)
                            conditions.energy.busses[bus.tag].power_draw   += conditions.energy.converters[converter.tag].inputs.power/bus.efficiency
                            conditions.energy.busses[bus.tag].current_draw  = conditions.energy.busses[bus.tag].power_draw/bus.voltage                            
                            
                        if isinstance(converter,RCAIDE.Library.Components.Powertrain.Converters.DC_Generator) or isinstance(converter,RCAIDE.Library.Components.Powertrain.Converters.PMSM_Generator):                              
                            compute_generator_performance(converter,conditions) 
                            conditions.energy.busses[bus.tag].power_draw   -= conditions.energy.converters[converter.tag].outputs.power/bus.efficiency
                            conditions.energy.busses[bus.tag].current_draw  = conditions.energy.busses[bus.tag].power_draw/bus.voltage                            
                        
        # ----------------------------------------------------------        
        # Section 3.0 Sources
        # ---------------------------------------------------------- 
        # 3.2 Fuel Sources  
        for fuel_line in fuel_lines:
            if fuel_line.active:
    
                # Update total mass flow of system   
                total_mdot  += conditions.energy.fuel_lines[fuel_line.tag].fuel_mass_flow_rate
                
                # Determine mass flow from each tank
                for tank in fuel_line.fuel_tanks:
                    tank.compute_tank_properties(state,fuel_line)  
                    state.conditions.energy.cumulative_fuel_consumption[1:,0] += np.cumsum(-np.diff(state.conditions.energy.fuel_lines.fuel_line.fuel_tanks[tank.tag].fuel_mass[:,0]))
                    
        # 3.2 Electric Sources 
        for bus in  busses:
            if bus.active: 
                # ------------------------------------------------------------------------------------------------------------------- 
                # 3.1 Batteries
                # -------------------------------------------------------------------------------------------------------------------       
                stored_results_flag       = False    
                for battery_module in  bus.battery_modules:                   
                    if bus.identical_battery_modules == False or stored_results_flag == False:
                        # run analysis  
                        stored_results_flag, stored_battery_cell_tag =  battery_module.energy_calc(state,bus,coolant_lines)
                    else:              
                        # use previous battery results 
                        battery_module.reuse_stored_data(state,bus,stored_results_flag, stored_battery_cell_tag)
                    
                # ------------------------------------------------------------------------------------------------------------------- 
                # 3.2 Fuel Cell Stacks
                # ------------------------------------------------------------------------------------------------------------------- 
                stored_results_flag       = False   
                stored_fuel_cell_tag      = None                  
                for fuel_cell_stack in  bus.fuel_cell_stacks:                   
                    if bus.identical_fuel_cell_stacks == False or stored_results_flag == False:
                        # run analysis  
                        stored_results_flag, stored_fuel_cell_tag =  fuel_cell_stack.energy_calc(state,bus,coolant_lines)
                    else:              
                        # use previous battery results 
                        fuel_cell_stack.reuse_stored_data(state,bus,stored_results_flag, stored_fuel_cell_tag)
                         
                    # compute mass flow rate                    
                    conditions.energy.busses[bus.tag].fuel_mass_flow_rate = state.conditions.energy.busses[bus.tag].fuel_cell_stacks[fuel_cell_stack.tag].H2_mass_flow_rate    
                          
                       
                    # Step 3: Compute bus properties          
                    bus.compute_distributor_conditions(state)
                    
                    # Step 4 : Battery Thermal Management Calculations                    
                    for coolant_line in coolant_lines:
                        for heat_exchanger in coolant_line.heat_exchangers: 
                            heat_exchanger.compute_heat_exchanger_performance(state,bus,coolant_line) 
                        for reservoir in coolant_line.reservoirs:   
                            reservoir.compute_reservior_coolant_temperature(state,coolant_line)
                                
                # Update total mass flow of system   
                total_mdot   += conditions.energy.busses[bus.tag].fuel_mass_flow_rate 
                                 
        if reverse_thrust ==  True:
            total_thrust =  total_thrust * -1    
            total_moment =  total_moment * -1                        
        conditions.energy.thrust_force_vector  = total_thrust
        conditions.energy.power                = total_mech_power 
        conditions.energy.thrust_moment_vector = total_moment 
        conditions.weights.vehicle_mass_rate   = total_mdot   
        
        return
    
    def unpack_unknowns(self,segment):
        """Unpacks the unknowns set in the mission to be available for the mission.
    
        Assumptions:
        N/A
        
        Source:
        N/A
        
        Inputs: 
            segment   - data structure of mission segment [-]
        
        Outputs: 
        
        Properties Used:
        N/A
        """            
         
        unknowns(segment)  
        for network in segment.analyses.energy.vehicle.networks:
            # Fuel unknowns 
            for fuel_line_i, fuel_line in enumerate(network.fuel_lines):    
                for propulsor_group in  fuel_line.assigned_propulsors:
                    propulsor = network.propulsors[propulsor_group[0]]
                    if propulsor.active: 
                        propulsor.unpack_unknowns(segment)
                #for fuel_tank in  fuel_line.fuel_tanks:
                    #fuel_tank.unpack_unknowns(fuel_line,segment)
                        
            # electric unknowns 
            for bus_i, bus in enumerate(network.busses):     
                for propulsor_group in  bus.assigned_propulsors:
                    propulsor = network.propulsors[propulsor_group[0]]
                    if propulsor.active: 
                        propulsor.unpack_unknowns(segment) 
        return    
     
    def residuals(self,segment): # these arenotusedin the mission solver per  seand needto berenamed
        """ This packs the residuals to be sent to the mission solver.
    
           Assumptions:
           None
    
           Source:
           N/A
    
           Inputs:
           state.conditions.energy:
               motor(s).torque                      [N-m]
               rotor(s).torque                      [N-m] 
           residuals soecific to the battery cell   
           
           Outputs:
           residuals specific to battery cell and network
    
           Properties Used: 
           N/A
       """         
        for network in segment.analyses.energy.vehicle.networks:
            for fuel_line_i, fuel_line in enumerate(network.fuel_lines):     
                for propulsor_group in  fuel_line.assigned_propulsors:
                    propulsor =  network.propulsors[propulsor_group[0]]
                    if propulsor.active: 
                        propulsor.pack_propulsor_residuals(segment) 
            for bus_i, bus in enumerate(network.busses):     
                for propulsor_group in  bus.assigned_propulsors:
                    propulsor =  network.propulsors[propulsor_group[0]]
                    if propulsor.active: 
                        propulsor.pack_propulsor_residuals(segment)   
        return      
     
# ----------------------------------------------------------------------
#  Component Container
# ---------------------------------------------------------------------- 
class Container(Component.Container):
    """ The Network container class 
    """
    def evaluate(self,state,center_of_gravity):
        """ This is used to evaluate the thrust and moments produced by the network.

            Assumptions:  
                If multiple networks are attached their performances will be summed

            Source:
                None 
        """ 
      
        self.evaluate(state,center_of_gravity)  
      


# ----------------------------------------------------------------------
#  Handle Linking
# ----------------------------------------------------------------------
Network.Container = Container