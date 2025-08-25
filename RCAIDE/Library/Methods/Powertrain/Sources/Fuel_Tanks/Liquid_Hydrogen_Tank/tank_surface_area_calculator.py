import numpy as np
from scipy.optimize import fsolve



# # Example usage
# r = 0.5     # radius
# l = 2       #****** l is the length of the cylinder
# h_true = 0.3
# v = l*(r**2*np.arccos((r - h_true)/r) - (r - h_true)*np.sqrt(2*r*h_true - h_true**2)) + np.pi*(2*r-h_true)**2 * (r-(2*r-h_true)/3)

# h_sol = solve_h_fsolve(r, l, v, h_guess=0.2)

# print("Target h =", h_true)
# print("Solved h =", h_sol)


# # Surface Area 

# area = np.pi*r*(h_sol) + 2*l*r*np.arccos((r-h_true)/r)
# print("Area=", area)
# #https://www.enground.com/Tank/Surface%20Area%20for%20Horizontal%20Cylindrical.html


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
    
