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
    return -0.2 * np.ones_like(m_l)

# ----------------------------
# Collocation residual
# ----------------------------
def residual(z, D_t, m_l0):
    cp = len(z)
    m_l = z
    dm_l = ode_mdot(m_l)  # or +0.2 if you want increasing

    R = D_t @ m_l - dm_l
    R[0] = m_l[0] - m_l0
    return R

def main():
    # -------------------------
    # Setup Chebyshev machinery
    # -------------------------
    segment  = RCAIDE.Framework.Mission.Segments.Segment()
    segment.state  = RCAIDE.Framework.Mission.Common.State()
    segment.state.numerics.number_of_control_points = 16
    differentials_dimensionless(segment)

    t0, tf = 0.0, 3600.0
    x = segment.state.numerics.dimensionless.control_points
    D = segment.state.numerics.dimensionless.differentiate
    D_t = D/(tf - t0)   # <-- make sure to include factor 2

    t_nodes = t0 + (tf - t0)*x   # physical times

    # -------------------------
    # Solve test IVP
    # -------------------------
    m_l0 = 1035.0
    guess = np.linspace(m_l0, m_l0 - 0.2*tf, len(x))  # linear guess

    sol = least_squares(residual, guess, args=(D_t, m_l0))
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