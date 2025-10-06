# RCAIDE/Methods/Thermal_Management/Batteries/Heat_Acquisition_System/Direct_Air_Heat_Acquisition/direct_convection_model.py
# 
# 
# Created:  Jul 2023, M. Clarke 

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  Compute Net Convected Heat 
# ---------------------------------------------------------------------------------------------------------------------- 
def air_cooled_performance(HAS,battery,_,coolant_line, Q_heat_gen,T_cell,T_scale,state):
    '''Computes the net heat removed by direct air heat acquisition system.

    Assumptions:
    1) Battery pack cell heat transfer can be modelled as a cooling columns in a cross-flow
    2) Isothermal battery cell - the temperature at the center of the cell is the same at 
    the surface of the cell

    Source: 
    Heat Transfer Model:
    Pakowski, Zdzisław. "Fundamentals of Heat and Mass Transfer, Frank P Incropera,
    David P DeWitt, Theodore L Bergman, Adrienne S Lavine, J. Wiley & Sons, Hoboken
    NJ (2007), 997 pp." (2007): 1683-1684.,  Chapter 7 pg 437-446 

    Inputs:  
              h                         (heat transfer coefficient)  [W/(m^2*K)] 
              As_cell                   (battery cell surface area)  [meters^2]
              H_cell                    (battery cell height)        [meters]
              T_ambient                 (ambient temperature)        [Kelvin]
              T_current                 (pack temperature)           [Kelvin]
              T_cell                    (battery cell temperature)   [Kelvin] 
              heat_transfer_efficiency                               [unitless]
              
      
      Outputs:
             T_cell                     (Updated battery cell temperature) [Kelvin]
  
    Properties Used:
    None 
    '''
    
    
    # battery properties  
    As_cell                  = battery.cell.surface_area 
    D_cell                   = battery.cell.diameter                     
    H_cell                   = battery.cell.height              
    cell_mass                = battery.cell.mass    
    Cp                       = battery.cell.specific_heat_capacity    
    Nn                       = battery.geometrtic_configuration.normal_count            
    Np                       = battery.geometrtic_configuration.parallel_count    
    n_total_module           = Nn*Np  
    h                        = HAS.convective_heat_transfer_coefficient 
    heat_transfer_efficiency = HAS.heat_transfer_efficiency   
    T_ambient                = state.conditions.freestream.temperature
       
    if n_total_module == 1: 
        # Using lumped model   
        Q_convec = h * As_cell * (T_cell - T_ambient)
        Q_heat_gen_tot = Q_heat_gen 

    else:   
        K_coolant   = state.conditions.freestream.thermal_conductivity
        nu_coolant  = state.conditions.freestream.kinematic_viscosity
        Pr_coolant  = state.conditions.freestream.prandtl_number
        rho_coolant = state.conditions.freestream.density
        Cp_coolant  = HAS.cooling_fluid.compute_cp(state.conditions.freestream.temperature,
                                                   state.conditions.freestream.pressure)
        V_coolant   = HAS.cooling_fluid.flowspeed  
        
        # Chapter 7 pg 437-446 of Fundamentals of heat and mass transfer 
        S_T = battery.geometrtic_configuration.normal_spacing          
        S_L = battery.geometrtic_configuration.parallel_spacing

        S_D = np.sqrt(S_T**2 + S_L**2)
        if 2*(S_D - D_cell) < (S_T - D_cell):
            V_max = V_coolant * (S_T / (2*(S_D - D_cell)))
        else:
            V_max = V_coolant * (S_T / (S_T - D_cell))

        T_film = (T_ambient + T_cell) / 2   
        Re_max = V_max * D_cell / nu_coolant   
        
        if np.all(Re_max > 10E2): 
            C = 0.35 * ((S_T/S_L)**0.2) 
            m = 0.6 
        else:
            C = 0.51
            m = 0.5  

        Pr_w_coolant = HAS.cooling_fluid.compute_prandtl_number(T_film)            
        Nu = C * (Re_max**m) * (Pr_coolant**0.36) * ((Pr_coolant/Pr_w_coolant)**0.25)           
        h = Nu * K_coolant / D_cell
        
        Tw_Ti = T_cell - T_ambient
        exp_term = np.exp((-np.pi * D_cell * n_total_module * h) / 
                         (rho_coolant * V_coolant * Nn * S_T * Cp_coolant))
        Tw_To = Tw_Ti * exp_term
        
        # Handle case where Tw_Ti is zero to avoid division by zero
        dT_lm = np.zeros_like(Tw_Ti)
        non_zero = Tw_Ti != 0.
        if np.any(non_zero):
            dT_lm[non_zero] = (Tw_Ti[non_zero] - Tw_To[non_zero]) / np.log(Tw_Ti[non_zero]/Tw_To[non_zero])
        
        Q_convec = heat_transfer_efficiency * h * np.pi * D_cell * H_cell * n_total_module * dT_lm 
        Q_heat_gen_tot = Q_heat_gen * n_total_module
    
    # Net heat and scaled temperature derivative
    Q_net = Q_heat_gen_tot - Q_convec  
    dT_dt = Q_net / (cell_mass * n_total_module * Cp)
    
    # Scale the temperature derivative
    dT_dt_scaled = dT_dt / T_scale

    state.conditions.energy.coolant_lines[coolant_line.tag][HAS.tag].total_heat_removed     = Q_convec*n_total_module

       
    
    return  dT_dt_scaled