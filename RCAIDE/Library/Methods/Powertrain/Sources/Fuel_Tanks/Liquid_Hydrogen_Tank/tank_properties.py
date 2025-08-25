import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

#from .tank_surface_area_calculator import compute_wetted_area

# ==================== GLOBAL  VARIABLES Constants (Will delete later) ====================
GRAV_CONST = 9.81
cv_g = 14300.0      # J/kg-K (gas)
cp_l = 9700.0       # J/kg-K (liquid)
P = 2e5             # Pa, assume fixed for now


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
    beta_l = 1.0 / max(T_liq, 1e-6)
    
    # Grashof-Prandtl product
    GrPr = (L_int**3 * rho_l**2 * GRAV_CONST * beta_l *
            abs(T_liq - T_int) * cp_l) / (mu_l * k_l + 1e-12)
    
    # Heat transfer coefficient
    alpha = C * (k_l / L_int) * (GrPr**n)
    
    # Heat flux
    Q = alpha * A_int * (T_liq - T_int)
    
    return Q

# ==================== Heat transfer between interface and gas ====================
def Q_gas_to_int(T_g, T_int, A_int, L_int, rho_g, cp_g, mu_g, k_g,
                 C=0.27, n=0.25):
    """
    Heat transfer from ullage to interface using Grashof–Prandtl natural convection.
    """

    beta_g = 1.0 / max(T_g, 1e-6)
    
    GrPr = (L_int**3 * rho_g**2 * GRAV_CONST * beta_g *
            abs(T_g - T_int) * cp_g) / (mu_g * k_g + 1e-12)
    
    alpha = C * (k_g / L_int) * (GrPr**n)
    
    Q = alpha * A_int * (T_g - T_int)
    
    return Q



# ==================== Tank ODEs ====================
def tank_odes(t, y,m_dot_g_out,m_dot_l_out): # Still need to bring in radius length of the tank 
    m_g, m_l, T_g, T_l, V_g, V_l = y
    
    # --- Interface saturation temperature ---
    T_int = PropsSI("T", "P", P, "Q", 0, "Hydrogen")  # [K]

    rho_l = PropsSI("D", "T", T_l, "Q", 0, "Hydrogen") # kg/m^3
    
    # --- Geometry placeholders ---
    A_int = 10.0   # m2
    L_int = 10   # m
    k_liq = PropsSI("L", "T", T_l, "Q", 0, "Hydrogen")   # thermal conductivity [W/m-K]
    mu_liq = PropsSI("V", "T", T_l, "Q", 0, "Hydrogen")  # viscosity [Pa·s]
    cp_liq = PropsSI("C", "T", T_l, "Q", 0, "Hydrogen")  # Cp [J/kg-K]

    # A_wet_liquid,A_wet_ullage = compute_wetted_area(r,l,V_g) # This should be volume of liquid 

    # --- Heat fluxes interface exchange ---
    Q_l_i = Q_liq_to_int(T_l, T_int, A_int, k_liq, L_int, 1.0, 0.25, cp_liq, mu_liq)
    Q_g_i = Q_gas_to_int(T_g, T_int, A_int, k_liq, L_int, 1.0, 0.25, cp_liq, mu_liq)

    # # --- Heat fluxes environment exchange ---
    # Q_e_g = Q_env_to_hydrogen(A_wet_ullage, T_h, T_g, N_layers=30)
    # Q_e_l = Q_env_to_hydrogen(A_wet_liquid, T_h, T_l, N_layers=30)

    Q_e_g = 20
    Q_e_l = 100
    
    # --- Enthalpies ---
    h_g = PropsSI("H", "T", T_int, "Q", 1, "Hydrogen")  # J/kg
    h_l = PropsSI("H", "T", T_int, "Q", 0, "Hydrogen")  # J/kg
    u_g = PropsSI("U", "T", T_int, "Q", 1, "Hydrogen")  # J/kg
    u_l = PropsSI("U", "T", T_int, "Q", 0, "Hydrogen")  # J/kg
    
    # --- Boil-off mass flow ---
    m_dot_bo = (Q_l_i + Q_g_i) / (h_g - h_l + 1e-9)
    
    # --- Mass balances ---
    dm_g = m_dot_bo - m_dot_g_out
    dm_l = -m_dot_bo - m_dot_l_out
    dV_g = -dm_l / rho_l
    dV_l = dm_l / rho_l
    
    # --- Energy balances ---
    dT_g = (-Q_g_i + Q_e_g - P*dV_g + dm_g*(h_g - u_g)) / (m_g*cv_g + 1e-9)
    dT_l = (-Q_l_i + Q_e_l- P*(-dV_g) + dm_l*(h_l - u_l)) / (m_l*cp_l + 1e-9)
    print(m_dot_bo)
    return [dm_g, dm_l, dT_g, dT_l, dV_g, dV_l]





# ==================== Main function ====================
def main():
    # Initial conditions
    y0 = [100, 400, 25.0, 20.0, 5,15]  # m_g, m_l, T_g, T_l, V_g,V_l
    t_span = (0, 3600)
    t_eval = np.linspace(*t_span, 500)

    # ==================== User input flows ====================
    m_dot_g_out = 0.01   # kg/s gas outflow
    m_dot_l_out = 0.1    # kg/s liquid outflow

    # Solve ODEs
    sol = solve_ivp(tank_odes, t_span, y0, t_eval=t_eval,args=(m_dot_g_out,m_dot_l_out))

    # --- Derived quantities ---
    T_ints = []
    m_boils = []
    Press = []
    for m_g, m_l, T_g, T_l, V_g,V_L in sol.y.T:
        P_now = m_g * 8.314/0.002016 * T_g / max(V_g,1e-9)  # ideal gas [Pa]
        # T_int = PropsSI("T", "P", P_now, "Q", 0, "Hydrogen")
        # h_g = PropsSI("H", "T", T_int, "Q", 1, "Hydrogen")
        # h_l = PropsSI("H", "T", T_int, "Q", 0, "Hydrogen")
        # m_dot_bo = (Q_l_i + Q_g_i) / (h_g - h_l + 1e-9)
        # T_ints.append(T_int)
        # #m_boils.append(m_dot_bo)
        # Press.append(P_now/1e5)  # bar

    # --- Plots ---
    plt.figure()
    plt.plot(sol.t, sol.y[2], label="Ullage Temp (T_g)")
    plt.plot(sol.t, sol.y[3], label="Liquid Temp (T_l)")
    #plt.plot(sol.t, T_ints, "--", label="Interface Temp (T_int)")
    plt.xlabel("Time [s]"); plt.ylabel("Temperature [K]")
    plt.legend(); plt.grid(True)

    # plt.figure()
    # plt.plot(sol.t, m_boils, label="Boil-off rate")
    # plt.xlabel("Time [s]"); plt.ylabel("m_dot_bo [kg/s]")
    # plt.legend(); plt.grid(True)

    plt.figure()
    #plt.plot(sol.t, Press, label="Tank Pressure")
    plt.xlabel("Time [s]"); plt.ylabel("Pressure [bar]")
    plt.legend(); plt.grid(True)

    plt.figure()
    plt.plot(sol.t, sol.y[0], label="Gas mass (m_g)")
    plt.plot(sol.t, sol.y[1], label="Liquid mass (m_l)")
    plt.xlabel("Time [s]"); plt.ylabel("Mass [kg]")
    plt.legend(); plt.grid(True)

    plt.show()

    
if __name__ == "__main__":
    main()
