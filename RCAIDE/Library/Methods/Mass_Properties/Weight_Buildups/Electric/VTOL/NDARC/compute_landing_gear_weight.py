# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/NDARC/compute_landing_gear_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  Compute landing gear weight
# ----------------------------------------------------------------------------------------------------------------------
# Direct port of Hydra's `src/Python/Stage_1/afdd/alighting.py::alighting_weight`
# (`github.com/VahanaOpenSource/vtol_sizing`). Unit-agnostic — a pure fraction of vehicle mass, no
# ft/lb conversion involved, unlike the wing methods.

_F_LG = 0.0497   # landing gear structure weight fraction of MTOW
_F_FA = 0.015    # landing gear fairing weight fraction of MTOW


def compute_landing_gear_weight(vehicle_mtow, tech_factor=1.0):
    """ Calculates the mass of the landing gear group (structure + fairing) as a fraction of
        vehicle MTOW.

        Source:
            Hydra `afdd/alighting.py::alighting_weight` —
            f_lg=0.0497 structure + f_fa=0.015 fairing, both x MTOW.

        Inputs:
            vehicle_mtow    vehicle max takeoff weight    [kg]
            tech_factor     technology weight-scaling factor, applied to the structural term only
                             (matches source — the fairing term is not tech-factor-scaled)    [Unitless]

        Outputs:
            weight:         dict with 'structure', 'fairing', 'total', all    [kg]
    """
    mass_lg = vehicle_mtow * _F_LG * tech_factor
    mass_fair = vehicle_mtow * _F_FA
    total = mass_lg + mass_fair
    return {'structure': mass_lg, 'fairing': mass_fair, 'total': total}
