# RCAIDE/Library/Missions/Common/Unpack_Unknowns/energy.py
# 
# 
# Created:  Jul 2023, M. Clarke
 
# ----------------------------------------------------------------------------------------------------------------------
#  Unpack Unknowns
# ---------------------------------------------------------------------------------------------------------------------- 
def unknowns(segment):  
    ACV_T      =  segment.assigned_control_variables.throttle
    ACV_TA     =  segment.assigned_control_variables.thrust_vector_angle
    ACV_RBPC   =  segment.assigned_control_variables.blade_pitch_command
    
    for network in segment.analyses.vehicle.networks: 
        if 'throttle' in segment: 
            for propulsor in network.propulsors: 
                segment.state.conditions.energy.propulsors[propulsor.tag].throttle[:,0] = segment.throttle
            
        if ACV_T.active:
            for i in range(len(ACV_T.assigned_propulsors)):
                propulsor_group = ACV_T.assigned_propulsors[i]
                for propulsor_name in propulsor_group:
                    # SPL fix (2026-07-23): was a bare reference assignment, which aliases this
                    # propulsor's condition.throttle array to the SAME object as
                    # state.unknowns["throttle_i"]. If any OTHER propulsor group in this segment
                    # is prescribed via segment.throttle (line 17's unconditional loop, which runs
                    # for every propulsor including this one), that in-place write then corrupts
                    # this group's solved unknown too, permanently, from the next solver iteration
                    # onward. .copy() breaks the aliasing.
                    segment.state.conditions.energy.propulsors[propulsor_name].throttle = segment.state.unknowns["throttle_" + str(i)].copy()
    
       # Thrust Vector Control 
        if ACV_TA.active:                
            for i in range(len(ACV_TA.assigned_propulsors)): 
                propulsor_group = ACV_TA.assigned_propulsors[i]
                for propulsor_name in propulsor_group:  
                    segment.state.conditions.energy.propulsors[propulsor_name].commanded_thrust_vector_angle = segment.state.unknowns["thrust_vector_angle_" + str(i)]
                     
        # Blade Pitch Command Control 
        if ACV_RBPC.active:                
            for i in range(len(ACV_RBPC.assigned_rotors)): 
                converter_group = ACV_RBPC.assigned_rotors[i] 
                for converter_name in converter_group:  
                    segment.state.conditions.energy.converters[converter_name].blade_pitch_command = segment.state.unknowns["blade_pitch_command_" + str(i)]                    
    return 
     
 
    
