import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from RCAIDE.Framework.Core import Units
import RCAIDE
from RCAIDE.Library.Mission.Common.Initialize import differentials_dimensionless
from RCAIDE.Library.Mission.Common.Initialize import time
from RCAIDE.Framework.Core import Data
from scipy.optimize import least_squares


# ----------------------------
# Simple test ODE
# ----------------------------
def ode_mdot(m_l):
    return np.array([[5.13380754],
       [5.13381948],
       [5.1339915 ],
       [5.13470254],
       [5.1364904 ],
       [5.13992852],
       [5.14548594],
       [5.15339711],
       [5.1635677 ],
       [5.17553605],
       [5.18849976],
       [5.20140396],
       [5.21307568],
       [5.22238114],
       [5.22837973],
       [5.23045137]])


# ----------------------------
# Collocation residual
# ----------------------------
def network_evaluate(z, D_t, m_l0):

    
    cp = len(z)
    m_l = z
    dm_l = ode_mdot(m_l)  # or +0.2 if you want increasing

    R = D_t @ m_l - dm_l[:,0]
    R[0] = m_l[0] - m_l0
    print(m_l)
    return R

def main():
    # -------------------------
    # Setup Chebyshev machinery
    # -------------------------
    segment  = RCAIDE.Framework.Mission.Segments.Segment()
    segment.state  = RCAIDE.Framework.Mission.Common.State()
    segment.state.numerics.number_of_control_points = 16
    differentials_dimensionless(segment)

    t0, tf = 0.0, 30
    x = segment.state.numerics.dimensionless.control_points
    D = segment.state.numerics.dimensionless.differentiate
    D_t = D/(tf - t0)   # <-- make sure to include factor 2

    t_nodes = t0 + (tf - t0)*x   # physical times

    # -------------------------
    # Solve test IVP
    # -------------------------
    m_l0 = 126607.
    guess = np.ones_like(x)*m_l0  # linear guess

    sol = least_squares(network_evaluate, guess[:,0], args=(D_t, m_l0))
    m_l = sol.x

    # -------------------------
    # Exact solution
    # -------------------------
    t_exact = np.linspace(t0, tf, 200)
    m_l_exact = m_l0 - 0.2*t_exact

    # -------------------------
    # Plot
    # -------------------------
    plt.plot(t_nodes, m_l, 'bo-', label="Collocation")
    plt.plot(t_exact, m_l_exact, 'r-', label="Exact")
    plt.xlabel("Time [s]")
    plt.ylabel("Liquid Mass [kg]")
    plt.title("Chebyshev IVP Test: Constant Outflow")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()