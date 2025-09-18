import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from RCAIDE.Framework.Core import Units
import RCAIDE
from RCAIDE.Library.Mission.Common.Initialize import differentials_dimensionless
from RCAIDE.Library.Mission.Common.Initialize import time
from RCAIDE.Framework.Core import Data
from scipy.optimize import least_squares
from scipy.optimize import fsolve
import scipy
import time
from RCAIDE.Framework.Optimization.Packages.scipy import scipy_setup





# ==================== Main function ====================
def main():
    method = 'least_squaresaaa'
    ti                   = time.time()

    # Initial Volume
    V_l = 29.379/2
    V_g = 1.6332/2
    # Inital Temperature
    T_l = 20
    T_g = 22

    # Inital Mass 
    rho_l = PropsSI("D", "T", T_l, "Q", 0, "Hydrogen") # kg/m^3
    rho_g = PropsSI("D", "T", T_g, "Q", 1, "Hydrogen") # kg/m^3

    m_g = rho_g * V_g
    m_l = rho_l* V_l



    # ==================== User input flows ====================


    # ==================== User input flows ====================
    m_dot_g_out = 0   # kg/s gas outflow
    m_dot_l_out = 0.2    # kg/s liquid outflow

    y0 = np.array((m_g, m_l, T_g, T_l , V_g, V_l))  # m_g, m_l, T_g, T_l, V_g,V_l




    segment  = RCAIDE.Framework.Mission.Segments.Segment()
    segment.state  = RCAIDE.Framework.Mission.Common.State()
    segment.state.numerics.number_of_control_points = 8
    differentials_dimensionless(segment)
    # segment.state.conditions.frames = Data()
    # segment.state.conditions.frames.inertial = Data()
    # segment.state.conditions.frames.inertial.time =  segment.state.ones_row(1)*0
    # segment.state.conditions.frames.inertial.time[:,0] = 

    t0, tf = 0.0, 3600.0

    x = segment.state.numerics.dimensionless.control_points
    D = segment.state.numerics.dimensionless.differentiate
    cp = segment.state.numerics.number_of_control_points
    y0 = build_initial_guess(m_g, m_l, T_g, T_l, V_g, V_l, tf, cp, m_dot_l_out)
    
    t_nodes = t0 + (tf - t0)*x
    D_t =D/(tf - t0)
    if method == 'least_squares':
        sol_cheb =  least_squares(residual, y0, args=(D_t,segment,m_dot_g_out,m_dot_l_out,y0), xtol=1e-6)
        cp = segment.state.numerics.number_of_control_points
        m_g = sol_cheb.x[:cp]
        m_l = sol_cheb.x[cp:2*cp]
        T_g = sol_cheb.x[2*cp:3*cp]
        T_l = sol_cheb.x[3*cp:4*cp]
        V_g = sol_cheb.x[4*cp:5*cp]
        V_l = sol_cheb.x[5*cp:6*cp]
    else:
        
        sol_cheb =  fsolve(residual, y0, args=(D_t,segment,m_dot_g_out,m_dot_l_out,y0), xtol=1e-6)
        cp = segment.state.numerics.number_of_control_points
        m_g = sol_cheb[:cp]
        m_l = sol_cheb[cp:2*cp]
        T_g = sol_cheb[2*cp:3*cp]
        T_l = sol_cheb[3*cp:4*cp]
        V_g = sol_cheb[4*cp:5*cp]
        V_l = sol_cheb[5*cp:6*cp]

    tf                   = time.time()
    elapsed_time         = ((tf-ti))
    print(' Simulation Time: ' + str(elapsed_time) + ' s')      

    # Plot T_g and T_l vs time
    plt.figure(figsize=(8,5))
    plt.plot(t_nodes, T_g, 'r-o', label="Gas Temperature $T_g$")
    plt.plot(t_nodes, T_l, 'b-s', label="Liquid Temperature $T_l$")

    plt.xlabel("Time [s]")
    plt.ylabel("Temperature [K]")
    plt.title("Gas and Liquid Temperature vs Time")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.figure(figsize=(8,5))
    # Gas mass
    plt.plot(t_nodes, m_g, 'g-o', label="Gas Mass $m_g$ [kg]")

    # Liquid mass

    plt.plot(t_nodes, m_l, 'b-s', label="Liquid Mass $m_l$ [kg]")

    # Labels and styling
    plt.xlabel("Time [s]")
    plt.ylabel("Mass [kg]")
    plt.title("Gas and Liquid Mass vs Time")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()



   


def residual(z,D_t,segment,m_dot_g_out,m_dot_l_out,y0):
    cp = segment.state.numerics.number_of_control_points
    
    m_g = z[:cp]
    m_l = z[cp:2*cp]
    T_g = z[2*cp:3*cp]
    T_l = z[3*cp:4*cp]
    V_g = z[4*cp:5*cp]
    V_l = z[5*cp:6*cp]


    m_g_0 = y0[:cp]
    m_l_0 = y0[cp:2*cp]
    T_g_0 = y0[2*cp:3*cp]
    T_l_0 = y0[3*cp:4*cp]
    V_g_0 = y0[4*cp:5*cp]
    V_l_0 = y0[5*cp:6*cp]

    dm_g, dm_l, dT_g, dT_l, dV_g, dV_l = tank_odes(m_g, m_l, T_g, T_l, V_g, V_l,m_dot_g_out,m_dot_l_out)

    R1 = D_t @ m_g - dm_g
    R1[0] = m_g[0] - m_g_0[0]

    R2 = D_t @ m_l - dm_l
    R2[0] = m_l[0] - m_l_0[0]

    
    R3 = D_t @ T_g - dT_g
    R3[0] = T_g[0] - T_g_0[0]
    
    R4 = D_t @ T_l - dT_l
    R4[0] = T_l[0] - T_l_0[0]
    
    R5 = D_t @ V_g - dV_g
    R5[0] = V_g[0] - V_g_0[0]
    
    R6 = D_t @ V_l - dV_l
    R6[0] = V_l[0] - V_l_0[0]

    return np.concatenate([R1, R2, R3, R4, R5, R6])

def build_initial_guess(m_g0, m_l0, T_g0, T_l0, V_g0, V_l0, tf, cp, m_dot_l_out=0.2):
    """
    Build a better collocation initial guess for the hydrogen tank.
    """
    # collocation nodes in [0,1]
    tau = np.linspace(0, 1, cp)

    # ---- Liquid mass (linearly decreasing) ----
    m_l_guess = m_l0 - m_dot_l_out * tf * tau

    # ---- Gas mass (roughly the boil-off replaces liquid outflow) ----
    m_g_guess = m_g0 + 0.05 * m_dot_l_out * tf * tau  # small growth

    # ---- Temperatures (small linear increase) ----
    T_g_guess = T_g0 + 1.0 * tau
    T_l_guess = T_l0 + 0.5 * tau

    # ---- Volumes (scale with masses/densities) ----
    rho_l0 = PropsSI("D", "T", T_l0, "Q", 0, "Hydrogen")
    rho_g0 = PropsSI("D", "T", T_g0, "Q", 1, "Hydrogen")

    V_l_guess = m_l_guess / rho_l0
    V_g_guess = m_g_guess / rho_g0

    # ---- Stack into vector [m_g, m_l, T_g, T_l, V_g, V_l] ----
    y0_guess = np.concatenate([
        m_g_guess, m_l_guess, T_g_guess, T_l_guess, V_g_guess, V_l_guess
    ])

    return y0_guess


# ==================== Tank ODEs ====================
def tank_odes( m_g, m_l, T_g, T_l, V_g, V_l,m_dot_g_out,m_dot_l_out): # Still need to bring in radius length of the tank or the whole fuel tank class
   

    # #Temperature of the interface assumes saturated hydrogen with same pressure as the ullage
    P = PropsSI("P", "T", T_g, "D", m_g / V_g, "Hydrogen")


    # --- Interface saturation temperature ---
    T_int = PropsSI("T", "P", P, "Q", 1, "Hydrogen")  # [K]
    print(T_int)

    rho_l = PropsSI("D", "T", T_l, "Q", 0, "Hydrogen") # kg/m^3
    
    # --- Geometry placeholders ---
    A_int = 16.06
    L_int = 1.542

    # Liquid properties
    k_liq = PropsSI("L", "T", T_l, "Q", 0, "Hydrogen")   # thermal conductivity [W/m-K]
    mu_liq = PropsSI("V", "T", T_l, "Q", 0, "Hydrogen")  # viscosity [Pa·s]
    cp_liq = PropsSI("C", "T", T_l, "Q", 0, "Hydrogen")  # Cp [J/kg-K]
    rho_l  = PropsSI("D", "T", T_l, "Q", 0, "Hydrogen")  # density [kg/m³]

    # Gas (ullage vapor) properties
    k_g   = PropsSI("L", "T", T_g, "Q", 1, "Hydrogen")   # thermal conductivity [W/m-K]
    mu_g  = PropsSI("V", "T", T_g, "Q", 1, "Hydrogen")   # viscosity [Pa·s]
    cp_g  = PropsSI("C", "T", T_g, "Q", 1, "Hydrogen")   # Cp [J/kg-K] # maybe wrong check thisn later
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
    # # --- Heat fluxes environment exchange ---
    # Q_e_g = Q_env_to_hydrogen(A_wet_ullage, T_h, T_g, N_layers=30)
    # Q_e_l = Q_env_to_hydrogen(A_wet_liquid, T_h, T_l, N_layers=30)

    Q_e_g = 20*np.ones_like(Q_l_i)
    Q_e_l = 25 *np.ones_like(Q_l_i)
    # --- Enthalpies ---
    h_g = PropsSI("H", "T", T_int, "Q", 1, "Hydrogen")  # J/kg
    h_l = PropsSI("H", "T", T_int, "Q", 0, "Hydrogen")  # J/kg
    u_g = PropsSI("U", "T", T_int, "Q", 1, "Hydrogen")  # J/kg
    u_l = PropsSI("U", "T", T_int, "Q", 0, "Hydrogen")  # J/kg
    
    # --- Natural Boil-off mass flow ---
    m_dot_bo = (Q_l_i + Q_g_i) / (h_g - h_l + 1e-9)

    # ---Vent(-) or boiled (+) mass flow rate ---
    m_dot_extra = m_g/(V_g*rho_l)*(m_dot_bo+m_dot_l_out) + m_dot_bo

    m_dot_extra_boi = np.where(m_dot_extra > 0, m_dot_extra, 0)
    m_dot_vent      = np.where(m_dot_extra < 0, -m_dot_extra, 0)

    m_dot_bo += m_dot_extra_boi
    m_dot_g_out += m_dot_vent

    # results["t"].append(t)
    # results["m_dot_bo_natural"].append(m_dot_bo-m_dot_extra_boi)
    # results["m_dot_additional"].append(m_dot_extra_boi)
    # results["m_dot_total"].append(m_dot_bo)
    # results["m_dot_vent"].append(m_dot_vent)

    # results["Pressure"].append(P)
    
    # --- Mass balances ---
    dm_g = m_dot_bo - m_dot_g_out
    dm_l = -m_dot_bo - m_dot_l_out
    dV_g = -dm_l / rho_l
    dV_l = dm_l / rho_l


    cv_g = 14300.0      # J/kg-K (gas)
    cp_l = 9700.0       # J/kg-K (liquid)
    
    # --- Energy balances ---
    dT_g = (-Q_g_i + Q_e_g - P*dV_g + dm_g*(h_g - u_g)) / (m_g*cv_g + 1e-9)
    dT_l = (-Q_l_i + Q_e_l- P*(dV_l) + dm_l*(h_l - u_l)) / (m_l*cp_l + 1e-9)

    return dm_g, dm_l, dT_g, dT_l, dV_g, dV_l


# ==================== Wetted Area Computer ====================
#****** l is the length of the cylinder
def equation(h, r, l, v):
    return l*(r**2 * np.arccos((r - h)/r) - (r - h)*np.sqrt(2*r*h - h**2)) + np.pi*(2*r-h)**2 * (r-(2*r-h)/3) - v

def solve_h_fsolve(r, t, v, h_guess=0.5):
    h_solution, = fsolve(equation, h_guess, args=(r, t, v))
    return h_solution

def compute_wetted_area(r,l,v):
    h_sol = solve_h_fsolve(r, l, v, h_guess=0.2)
    total_area = 2*np.pi*r*l + 4*np.pi*r**2
    area_liquid = np.pi*r*(h_sol) + 2*l*r*np.arccos((r-h_sol)/r)
    area_ullage = total_area-area_liquid
    return area_liquid, area_ullage
    


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
    GrPr = (L_int**3 * rho_l**2 * 9.81 * beta_l *
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

    beta_g = 1.0 / T_g
    
    GrPr = (L_int**3 * rho_g**2 * 9.81 * beta_g *
            abs(T_g - T_int) * cp_g) / (mu_g * k_g + 1e-12)
    
    alpha = C * (k_g / L_int) * (GrPr**n)
    
    Q = alpha * A_int * (T_g - T_int)
    
    return Q


if __name__ == "__main__":
    main()
    plt.show()