import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

# ---------------------------------------------------
# Parameters
# ---------------------------------------------------
dt = 10.0         # time step [s]
t_end = 10*3600    # total time [s] = 3 hrs
n_steps = int(t_end/dt)

fluid = "Hydrogen"

# Tank geometry
V_total = 30.0        # m^3 total tank volume
fill_fraction = 0.5
rhoL0 = CP.PropsSI("D","T",20,"Q",0,fluid)   # liquid density at 20K
rhoV0 = CP.PropsSI("D","T",20,"Q",1,fluid)   # vapor density at 20K

m_liquid = fill_fraction * V_total * rhoL0
m_vapor  = (1-fill_fraction) * V_total * rhoV0

# Heat leaks
Q_liquid = 50.0   # W
Q_vapor  = 10.0   # W

# Outflows
mdot_fuel = 0.5      # fuel outflow [kg/s]

# Target pressure (constant venting regime)
P_set = 1.0e5   # Pa (1 bar)

# ---------------------------------------------------
# Storage arrays
# ---------------------------------------------------
time_hist, P_hist, TL_hist, TV_hist = [], [], [], []
VL_hist, VV_hist = [], []

# ---------------------------------------------------
# Time stepping
# ---------------------------------------------------
for step in range(n_steps):
    t = step*dt
    
    # Current liquid and vapor volumes
    V_liquid = m_liquid / rhoL0
    V_vapor  = V_total - V_liquid
    
    # Pressure fixed at setpoint
    P = P_set
    
    # Saturation temp at P_set
    T_sat = CP.PropsSI("T","P",P,"Q",0,fluid)
    
    # Liquid & vapor temps
    TL = T_sat
    TV = T_sat + 2.0   # assume ullage superheated by 2 K
    
    # Latent heat
    h_fg = CP.PropsSI("H","P",P,"Q",1,fluid) - CP.PropsSI("H","P",P,"Q",0,fluid)
    
    # Evaporation from heat leak
    Q_interface = Q_liquid + Q_vapor
    mdot_evap = Q_interface / h_fg
    
    # --- Mass balances ---
    dm_liquid = -mdot_evap * dt
    mdot_vent = max(mdot_evap - mdot_fuel, 0.0)  # venting to keep pressure constant
    dm_vapor = (mdot_evap - mdot_vent - mdot_fuel) * dt
    
    # update masses
    m_liquid = max(m_liquid + dm_liquid, 0.0)
    m_vapor  = max(m_vapor + dm_vapor, 0.0)
    
    # recompute volumes
    V_liquid = m_liquid / rhoL0
    V_vapor  = V_total - V_liquid
    
    # store
    time_hist.append(t/3600) # hr
    P_hist.append(P/1e5)     # bar
    TL_hist.append(TL)
    TV_hist.append(TV)
    VL_hist.append(V_liquid)
    VV_hist.append(V_vapor)

# ---------------------------------------------------
# Plotting
# ---------------------------------------------------
plt.figure(figsize=(10,8))

plt.subplot(3,1,1)
plt.plot(time_hist,P_hist)
plt.ylabel("Pressure [bar]")

plt.subplot(3,1,2)
plt.plot(time_hist,TL_hist,label="Liquid T")
plt.plot(time_hist,TV_hist,label="Ullage T")
plt.ylabel("Temperature [K]")
plt.legend()

plt.subplot(3,1,3)
plt.plot(time_hist,VL_hist,label="Liquid volume")
plt.plot(time_hist,VV_hist,label="Ullage volume")
plt.ylabel("Volume [m³]")
plt.xlabel("Time [hr]")
plt.legend()

plt.tight_layout()
plt.show()
