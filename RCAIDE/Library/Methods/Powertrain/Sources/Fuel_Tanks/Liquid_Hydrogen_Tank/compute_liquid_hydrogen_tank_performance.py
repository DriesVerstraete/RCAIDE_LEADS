# RCAIDE/Methods/Powertrain/Sources/Fuel_Tanks/
# 
# 
# Created:  
# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports

# Python imports
from copy import deepcopy
import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.optimize import brentq
# ----------------------------------------------------------------------------------------------------------------------
#  METHOD
# ----------------------------------------------------------------------------------------------------------------------  
def compute_liquid_hydrogen_tank_performance(fuel_tank,state,distributor):

    distributor_conditions = state.conditions.energy.fuel_lines[distributor.tag]          
    tank_conditions       = distributor_conditions.fuel_tanks[fuel_tank.tag]   
    
    m_g = state.unknowns.network[fuel_tank.tag + '_ullage_mass']
    m_l = state.unknowns.network[fuel_tank.tag + '_liquid_mass']
    T_g = state.unknowns.network[fuel_tank.tag + '_ullage_temperature']
    T_l = state.unknowns.network[fuel_tank.tag + '_liquid_temperature']
    V_g = state.unknowns.network[fuel_tank.tag + '_ullage_volume']
    V_l = state.unknowns.network[fuel_tank.tag + '_liquid_volume'] 

    fuel_flow_split_ratio  = fuel_tank.fuel_selector_valve.fuel_flow_split_ratio 
    if fuel_tank.xz_plane_symmetric:
        fuel_flow_split_ratio  /=2

    tank_conditions.mass_flow_rate  =  distributor_conditions.fuel_mass_flow_rate * fuel_flow_split_ratio    
    
    m_dot_l_out =  tank_conditions.mass_flow_rate  
    m_dot_g_out =  tank_conditions.vent_rate

    # --- Interface saturation Pressure ---
    P = PropsSI("P", "T", T_g[:,0], "D", m_g[:,0] / V_g[:,0], "Hydrogen")


    # --- Interface saturation temperature ---
    T_int = PropsSI("T", "P", P, "Q", 1, "Hydrogen")  # [K]  
    
    # --- Geometry placeholders ---
    L_int,A_int= compute_interface_geometric_properties(fuel_tank, V_l[:,0]) 

    # Liquid properties
    k_liq = PropsSI("L", "T", T_l, "Q", 0, "Hydrogen")   # thermal conductivity [W/m-K]
    mu_liq = PropsSI("V", "T", T_l, "Q", 0, "Hydrogen")  # viscosity [Pa·s]
    cp_liq = PropsSI("C", "T", T_l, "Q", 0, "Hydrogen")  # Cp [J/kg-K]
    rho_l  = PropsSI("D", "T", T_l, "Q", 0, "Hydrogen")  # density [kg/m³]

    # Gas (ullage vapor) properties
    k_g   = PropsSI("L", "T", T_g, "Q", 1, "Hydrogen")   # thermal conductivity [W/m-K]
    mu_g  = PropsSI("V", "T", T_g, "Q", 1, "Hydrogen")   # viscosity [Pa·s]
    cp_g  = PropsSI("Cpmass", "T", T_g, "Q", 1, "Hydrogen")   # Cp [J/kg-K] # maybe wrong check thisn later
    rho_g = PropsSI("D", "T", T_g, "Q", 1, "Hydrogen")   # density [kg/m³]

    # --- Heat fluxes interface exchange ---
    Q_l_i = Q_liq_to_int(
        T_l, T_int, A_int, L_int,
        rho_l, cp_liq, mu_liq, k_liq
    )
    
    Q_g_i = Q_gas_to_int(
        T_g, T_int, A_int, L_int,
        rho_g, cp_g, mu_g, k_g
    )

    # --- Heat fluxes environment exchange ---
    # Q_e_g = Q_env_to_hydrogen(A_wet_ullage, T_h, T_g, N_layers=30)
    # Q_e_l = Q_env_to_hydrogen(A_wet_liquid, T_h, T_l, N_layers=30)
    
    # --- Heat fluxes environment exchange assume constant for testing code ---
    Q_e_g = 200*np.ones_like(Q_l_i)
    Q_e_l = 200*np.ones_like(Q_l_i)

    # --- Enthalpies ---
    h_g = PropsSI("H", "T", T_int, "Q", 1, "Hydrogen")  # J/kg
    h_l = PropsSI("H", "T", T_int, "Q", 0, "Hydrogen")  # J/kg
    u_g = PropsSI("U", "T", T_int, "Q", 1, "Hydrogen")  # J/kg
    u_l = PropsSI("U", "T", T_int, "Q", 0, "Hydrogen")  # J/kg
    
    # --- Natural Boil-off mass flow ---
    m_dot_bo = (Q_l_i + Q_g_i) / (h_g - h_l + 1e-9)

    # ---Vent(-) or boiled (+) mass flow rate ---
    m_dot_extra = m_g[:,0]/(V_g[:,0]*rho_l)*(m_dot_bo+ m_dot_l_out[:,0]) + m_dot_bo

    m_dot_extra_boi = np.where(m_dot_extra > 0, m_dot_extra, 0)
    m_dot_vent      = np.where(m_dot_extra < 0, -m_dot_extra, 0)

    m_dot_bo         += m_dot_extra_boi
    m_dot_g_out[:,0] += m_dot_vent
    
    # --- Mass balances ---
    dm_g = m_dot_bo - m_dot_g_out[:,0]
    dm_l = -m_dot_bo - m_dot_l_out[:,0]
    dV_g = -dm_l / rho_l
    dV_l = dm_l / rho_l
    
    # --- Energy balances ---
    dT_g = (-Q_g_i + Q_e_g - P*dV_g + dm_g*(h_g - u_g)) / (m_g[:,0]*cp_g + 1e-9)
    dT_l = (-Q_l_i + Q_e_l- P*(dV_l) + dm_l*(h_l - u_l)) / (m_l[:,0]*cp_liq + 1e-9)
 
    D = state.numerics.time.differentiate

    state.residuals.network[fuel_tank.tag + '_ullage_mass']           = np.dot(D, m_g)[:, 0] - dm_g[:, 0]
    state.residuals.network[fuel_tank.tag + '_ullage_mass'][0]        = m_g[0] -  tank_conditions.ullage_mass[0] 
    state.residuals.network[fuel_tank.tag + '_liquid_mass']           = np.dot(D, m_l)[:, 0] - dm_l[:, 0]
    state.residuals.network[fuel_tank.tag + '_liquid_mass'][0]        = m_l[0] - tank_conditions.mass[0] 
    state.residuals.network[fuel_tank.tag + '_ullage_temperature']    = np.dot(D, T_g)[:, 0] - dT_g[:, 0]
    state.residuals.network[fuel_tank.tag + '_ullage_temperature'][0] = T_g[0] - tank_conditions.ullage_temperature[0,0] 
    state.residuals.network[fuel_tank.tag + '_liquid_temperature']    = np.dot(D, T_l)[:, 0] - dT_l[:, 0]
    state.residuals.network[fuel_tank.tag + '_liquid_temperature'][0] = T_l[0] - tank_conditions.liquid_temperature[0,0] 
    state.residuals.network[fuel_tank.tag + '_ullage_volume']         = np.dot(D, V_g)[:, 0] - dV_g[:, 0]
    state.residuals.network[fuel_tank.tag + '_ullage_volume'][0]      = V_g[0] - tank_conditions.ullage_volume[0,0] 
    state.residuals.network[fuel_tank.tag + '_liquid_volume']         = np.dot(D, V_l)[:, 0] - dV_l[:, 0]
    state.residuals.network[fuel_tank.tag + '_liquid_volume'][0]      = V_l[0] - tank_conditions.liquid_volume[0,0] 

    tank_conditions.ullage_mass[1:,0]        = state.unknowns.network[fuel_tank.tag + '_ullage_mass'][1:,0]
    tank_conditions.mass[1:,0]               = state.unknowns.network[fuel_tank.tag + '_liquid_mass'][1:,0]
    tank_conditions.ullage_temperature[1:,0] = state.unknowns.network[fuel_tank.tag + '_ullage_temperature'][1:,0]
    tank_conditions.liquid_temperature[1:,0] = state.unknowns.network[fuel_tank.tag + '_liquid_temperature'][1:,0]
    tank_conditions.ullage_volume[1:,0]      = state.unknowns.network[fuel_tank.tag + '_ullage_volume'][1:,0]
    tank_conditions.liquid_volume[1:,0]      = state.unknowns.network[fuel_tank.tag + '_liquid_volume'][1:,0]
    tank_conditions.vent_rate                = m_dot_g_out
    tank_conditions.boil_off_rate[:,0]       = m_dot_bo
    
    if fuel_tank.xz_plane_symmetric:
        symmetric_tag = fuel_tank.tag  + "_symmetric"
        distributor_conditions.fuel_tanks[symmetric_tag] = deepcopy(distributor_conditions.fuel_tanks[fuel_tank.tag])

    return

# =========================================================================================    

# ==================== Heat transfer between Hydrogen and environment ====================
def Q_env_to_hydrogen(A_wet, T_h, T_c, N_layers=30):

    # Default constants from Keller et al.
    C_r = 5.39e-10             # Radiation constant
    emittance = 0.031          # Effective emittance
    C_s = 8.95e-8              # Solid conduction constant
    layer_density = 30 * 100   # 30 layers/cm → 3000 layers/m
    C_g = 1.46e4               # Gas conduction constant
    # Convert 1 torr = 133.322 Pa
    P = 1e-6 * 133.322         # Vacuum pressure [Pa]

    # Heat flux due to radiation
    q_rad = C_r * emittance / N_layers * (T_h**4.67 - T_c**4.67)

    # Heat flux due to solid conduction
    q_solid_cond = (C_s * layer_density**2.56 / N_layers *
                    (T_h + T_c) / 2 * (T_h - T_c))

    # Heat flux due to gas conduction
    q_gas_cond = C_g * P / N_layers * (T_h**0.52 - T_c**0.52)

    # Total MLI heat flux
    q_MLI = q_rad + q_solid_cond + q_gas_cond

    # Total heat transfer rate
    Q_dot = A_wet * q_MLI

    return Q_dot

# ==================== Heat transfer between interface and liquid ====================
def Q_liq_to_int(T_liq, T_int, A_int, L_int, rho_l, cp_l, mu_l, k_l,
                 C=0.27, n=0.25):
    """
    Heat transfer from liquid to interface using Grashof–Prandtl natural convection.
    """
    # Thermal expansion coefficient (approx. ideal fluid ~ 1/T)
    beta_l = 1.0 / T_liq
    
    # Grashof-Prandtl product
    GrPr = (L_int**3 * rho_l**2 * 9.81 * beta_l[:,0] *
            abs(T_liq[:,0] - T_int) * cp_l) / (mu_l * k_l + 1e-12)
    
    # Heat transfer coefficient
    alpha = C * (k_l / L_int) * (GrPr**n)
    
    # Heat flux
    Q = alpha * A_int * (T_liq[:,0] - T_int)
    
    return Q

# ==================== Heat transfer between interface and gas ====================
def Q_gas_to_int(T_g, T_int, A_int, L_int, rho_g, cp_g, mu_g, k_g,
                 C=0.27, n=0.25):
    """
    Heat transfer from ullage to interface using Grashof Prandtl natural convection.
    """

    beta_g = 1.0 / T_g
    
    GrPr = (L_int**3 * rho_g**2 * 9.81 * beta_g[:,0] *
            abs(T_g[:,0] - T_int) * cp_g) / (mu_g * k_g + 1e-12)
    
    alpha = C * (k_g / L_int) * (GrPr**n)
    
    Q = alpha * A_int * (T_g[:,0] - T_int)
    
    return Q

# ==================== Compute Height of the Liquid in the tank given a particular volume ====================

def solve_height_equation(h, r, l, v):
    """
    v = (
    l * (
        r**2 * np.arccos((r - h) / r)
        - (r - h) * np.sqrt(2 * r * h - h**2)
    )
    + (np.pi / 3) * h**2 * (3 * r - h)
)
    """

    return  l * (r**2 * np.arccos((r - h) / r) - (r - h) * np.sqrt(2 * r * h - h**2))+ (np.pi / 3) * h**2 * (3 * r - h) - v


def compute_liquid_height_per_volume(r, l, v):
    """Scalar: solve for liquid height given a single volume v."""
    return brentq(solve_height_equation, 0.0, 2*r, args=(r, l, v))

def compute_liquid_height(r, l, v_array):
    """Vectorized: works for scalar or array input v_array."""
    v_array = np.atleast_1d(v_array)
    h_array = np.fromiter(
        (compute_liquid_height_per_volume(r, l, vv) for vv in v_array),
        dtype=float
    )
    return h_array.reshape(np.shape(v_array))

def compute_interface_geometric_properties(fuel_tank, v_l):
    l_in = fuel_tank.inner_length
    r_in = fuel_tank.inner_diameter / 2
    v_l[v_l<0] = 0

    # h is array if v_l is array
    h = compute_liquid_height(r_in, l_in, v_l)

    lamda = h / (2 * r_in)
    L_int = 4 * r_in * np.sqrt(lamda - lamda**2)
    A_int = (np.pi * 2 * r_in * L_int**2 / 4) + L_int * l_in
    return L_int, A_int