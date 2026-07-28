# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Conventional/General_Aviation/Roskam/__init__.py
#

"""RCAIDE Package Setup

Roskam, Airplane Design Parts I Through VIII (2nd ed., DARcorp., 2003) GA weight-estimation
relationships — ported 2026-07-28 from Chakraborty & Mishra, Journal of Aircraft, Vol. 58, No. 4,
2021, Appendix C, Table C1 (confirmed by direct visual read of the paper's equation image; not
independently cross-checked against Roskam's own textbook, which is not available locally, unlike
Raymer's — see 20-rcaide-weight-method-porting-inventory.md for full detail).

Sibling of `Raymer`/`FLOPS` — not a standalone selectable `weights_analysis.method` on its own;
exists to supply the Roskam side of Chakraborty2022's per-component averaging (each component
where the paper cites both Raymer and Roskam sources).
"""

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

from .compute_fuselage_weight       import compute_fuselage_weight
from .compute_wing_weight           import compute_wing_weight_1, compute_wing_weight_2
from .compute_vertical_tail_weight  import compute_vertical_tail_weight
from .compute_electrical_weight     import compute_electrical_weight
from .compute_furnishings_weight    import compute_furnishings_weight_1, compute_furnishings_weight_2
