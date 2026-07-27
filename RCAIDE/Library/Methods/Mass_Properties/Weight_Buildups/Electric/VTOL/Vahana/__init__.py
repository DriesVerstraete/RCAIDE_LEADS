# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/Vahana/__init__.py
#

"""RCAIDE Package Setup

`Vahana` is an alias for `Physics_Based` (Project Vahana conceptual trade study method) — every
name here re-exports from `Physics_Based` directly so the two methods can never drift apart.
"""

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

from ..Physics_Based.compute_boom_weight            import compute_boom_weight
from ..Physics_Based.compute_fuselage_weight        import compute_fuselage_weight
from .compute_operating_empty_weight                import compute_operating_empty_weight
from ..Physics_Based.compute_wing_weight            import compute_wing_weight
from ..Physics_Based.dynamo_supply_mass_estimation  import dynamo_supply_mass_estimation
