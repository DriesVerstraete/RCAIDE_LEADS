# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/Chakraborty/__init__.py
#

"""RCAIDE Package Setup

Airframe-structure weight method ported from Chakraborty, I., and Mishra, A. A., "Generalized
Energy-Based Flight Vehicle Sizing and Performance Analysis Methodology," Journal of Aircraft,
Vol. 58, No. 4, 2021 — General Aviation statistical weight-estimation relationships (Raymer +
Roskam), averaged per component where the paper cites both sources, with a tilt-mechanism weight
penalty (`KTW`, default 1.20) applied to wing/horizontal-tail/canard.

**Not physics-based** — statistical regressions fit to existing GA aircraft, included as an
independent cross-check, not because this method is expected to be more physically correct than
`Vahana`/`NDARC`/`Hydra`.

Covers wing, horizontal tail, vertical tail, fuselage, landing gear, electrical, and furnishings
only — everything else (motor, rotor, BRS, battery, seats/avionics/ECS, wiring, thermal
management) falls back to `Vahana`/`Physics_Based`, same role it plays for `NDARC`/`Hydra`/
`CADDEE_LPC`'s own gaps. Fuel system/flight controls/hydraulics come from Raymer's own bundled
`compute_systems_weight` (all three now source-corrected — see the inventory doc).

See `01-mission-profiles/01-docs/rcaide/20-rcaide-weight-method-porting-inventory.md` for the
full Table C1 breakdown, the Raymer-textbook cross-check, and the RCAIDE bugs found/fixed along
the way.

See Also
--------
RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Conventional.General_Aviation.Raymer
RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Conventional.General_Aviation.Roskam
RCAIDE.Library.Methods.Mass_Properties.Weight_Buildups.Electric.VTOL.Vahana
"""

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

from .compute_operating_empty_weight import compute_operating_empty_weight
