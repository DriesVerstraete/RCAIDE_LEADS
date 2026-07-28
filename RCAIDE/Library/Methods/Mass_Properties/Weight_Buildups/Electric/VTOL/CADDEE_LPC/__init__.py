# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/CADDEE_LPC/__init__.py
#

"""RCAIDE Package Setup

Airframe-structure weight method ported from LSDOlab/UCSD's `CADDEE_alpha`
(`github.com/LSDOlab/CADDEE_alpha`, `core/aircraft/models/weights/nasa_lpc/m4_nasa_lpc.py`) — a
linear-regression surrogate fit against a DOE sweep of the NASA Lift+Cruise (LPC) reference
vehicle's `aframe`-FEA-optimized structure (Ruh & Gandhi, AIAA 2023-4280).

**Valid ONLY for lift+cruise-class vehicles**, and only near the LPC reference vehicle's own
design-space neighborhood — this is NOT a general VTOL weight method the way `Vahana`/`NDARC`/
`Hydra` are. Do not select `method='CADDEE_LPC'` for tiltrotor, tiltwing, or tilt-stopped-rotor
vehicles.

Covers wing, fuselage, boom, and empennage mass only — CADDEE_alpha has no motor, rotor, BRS, or
battery weight model at all; everything else falls back to a single decoupled `Vahana` pass (see
`compute_operating_empty_weight.py`'s module docstring for why it's decoupled rather than re-run
every iteration).

See `01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md` for the full
validation record and numeric cross-check against `Vahana`/`NDARC`/`Hydra`/real-NDARC-Fortran on
the real `lift_cruise` vehicle.

See Also
--------
RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.NDARC
RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.Hydra
RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.Vahana
"""

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

from .compute_wing_weight import compute_wing_weight
from .compute_fuselage_weight import compute_fuselage_weight
from .compute_boom_weight import compute_boom_weight
from .compute_empennage_weight import compute_empennage_weight
from .compute_operating_empty_weight import compute_operating_empty_weight
