# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/Hydra/__init__.py
#

"""RCAIDE Package Setup

Collection of aircraft weight estimation methods based on direct ports of Hydra's class-based
weight layer (`src/Python/Stage_1/wing_class.py`, `blade_wt_modelv2.py`, `fuselage_v2.py` and
supporting geometry helpers, `github.com/VahanaOpenSource/vtol_sizing`) — the layer SUAVE/RCAIDE's
`empty_hydra.py` (audited 2026-07-15) was originally ported from. Generally simpler and less
tiltrotor-differentiated than the `NDARC` method. See
`01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md` for the full
component-by-component porting inventory and status.

See Also
--------
RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.NDARC
RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.Vahana
RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.Physics_Based
"""

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

from .compute_wing_weight                        import compute_wing_weight_group
from .compute_motor_weight                        import compute_motor_weight
from .compute_rotor_weight                        import compute_rotor_weight
from .compute_fuselage_weight                     import compute_fuselage_weight
from .compute_operating_empty_weight              import compute_operating_empty_weight

# BRS and landing gear: Hydra's own Stage_1 layer does not reimplement these separately, it calls
# the same `afdd/emergency_sys.py` / `afdd/alighting.py` formulas the `NDARC` method already ports
# bit-for-bit. Re-exported directly (not duplicated) so the two can never drift apart, following
# the same alias pattern as `Vahana`/`Physics_Based`.
from ..NDARC.compute_emergency_system_weight     import compute_emergency_system_weight
from ..NDARC.compute_landing_gear_weight         import compute_landing_gear_weight
