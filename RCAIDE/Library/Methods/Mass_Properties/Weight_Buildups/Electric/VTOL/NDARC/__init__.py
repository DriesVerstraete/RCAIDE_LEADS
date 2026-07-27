# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/NDARC/__init__.py
#

"""RCAIDE Package Setup

Collection of aircraft weight estimation methods based on direct ports of NDARC's published
AFDD82/83/84/00 statistical regressions (Hydra's `src/Python/Stage_1/afdd/*.py`,
`github.com/VahanaOpenSource/vtol_sizing`). See
`01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md` for the full
component-by-component porting inventory and status.

See Also
--------
RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.Hydra
RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.Vahana
RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.Physics_Based
"""

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

from .compute_wing_weight import compute_wing_weight_group, compute_wing_tip_mass
from .compute_emergency_system_weight import compute_emergency_system_weight
from .compute_landing_gear_weight import compute_landing_gear_weight
from .compute_flight_control_system_weight import compute_flight_control_system_weight
