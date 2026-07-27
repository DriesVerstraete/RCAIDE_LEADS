# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/NDARC/compute_emergency_system_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  Compute emergency system (BRS) weight
# ----------------------------------------------------------------------------------------------------------------------
# Direct port of Hydra's `src/Python/Stage_1/afdd/emergency_sys.py::emergency_sys`
# (`github.com/VahanaOpenSource/vtol_sizing`). Unit-agnostic — a pure fraction of vehicle mass, no
# ft/lb conversion involved, unlike the wing methods.

def compute_emergency_system_weight(vehicle_mtow, tech_factor=1.0):
    """ Calculates the mass of a ballistic recovery system (BRS) as a fraction of vehicle MTOW.

        Source:
            Hydra `afdd/emergency_sys.py::emergency_sys` — 0.0233x MTOW.

        Inputs:
            vehicle_mtow    vehicle max takeoff weight    [kg]
            tech_factor     technology weight-scaling factor    [Unitless]

        Outputs:
            mass_brs:       ballistic recovery system mass    [kg]
    """
    mass_brs = 0.0233 * vehicle_mtow * tech_factor
    return mass_brs
