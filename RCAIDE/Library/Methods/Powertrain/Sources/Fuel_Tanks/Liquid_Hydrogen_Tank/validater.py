import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

# -------------------------------------------------------------------
# Parameters for Newton’s cooling
# -------------------------------------------------------------------
h   = 100.0      # W/m^2-K
A   = 0.1       # m^2
m   = 2.0       # kg
cp  = 900.0     # J/kg-K
T_inf = 25.0    # °C
k = h*A/(m*cp)  # 1/s

T0 = 100.0      # initial temperature
t0, tf = 0.0, 3600.0

# -------------------------------------------------------------------
# ODE solver version
# -------------------------------------------------------------------
def cooling_ode(t, T):
    return -(h*A/(m*cp))*(T - T_inf)

t_eval = np.linspace(t0, tf, 200)
sol_ivp = solve_ivp(cooling_ode, (t0, tf), [T0], t_eval=t_eval)

# -------------------------------------------------------------------
# Chebyshev collocation version
# -------------------------------------------------------------------
def chebyshev_data(N=8, integration=True):
    N = int(N)
    if N <= 1:
        raise ValueError("N must be >= 2")
    x = 0.5*(1 - np.cos(np.pi*np.arange(0, N)/(N-1)))  # nodes in [0,1]
    c = np.ones(N); c[0] = c[-1] = 2.0
    c = c * ((-1.0)**np.arange(N))
    X  = np.tile(x, (N,1)).T
    dX = X - X.T + np.eye(N)
    D = (np.outer(c,1/c)) / dX
    D = D - np.diag(np.sum(D.T, axis=0))
    if not integration:
        return x, D, None
    I_int = np.linalg.inv(D[1:,1:])
    I = np.vstack([np.zeros((1, N-1)), I_int])
    I = np.hstack([np.zeros((N,1)), I])
    return x, D, I

N = 32
x, D, I = chebyshev_data(N)
t_nodes = t0 + (tf - t0)*x
D_t = D/(tf - t0)

def rhs(T):
    return -k*(T - T_inf)

def residual(z):
    T = z
    R = D_t @ T - rhs(T)
    R[0] = T[0] - T0
    return R

z0 = np.full(N, T0)
sol_cheb = least_squares(residual, z0, method="trf", xtol=1e-12, ftol=1e-12)
T_cheb = sol_cheb.x

# -------------------------------------------------------------------
# Plot comparison
# -------------------------------------------------------------------
plt.figure(figsize=(7,5))
plt.plot(sol_ivp.t, sol_ivp.y[0], label="solve_ivp (RK45)", lw=2)
plt.plot(t_nodes, T_cheb, "o-", label="Chebyshev collocation", lw=1.5)
plt.axhline(T_inf, color="k", ls="--", label="Ambient")
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.title("Cooling of a Hot Block: solve_ivp vs. Chebyshev")
plt.legend(); plt.grid(True)
plt.show()
