# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/Vahana/compute_operating_empty_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

# `Vahana` is an alias for `Physics_Based` — same Project Vahana conceptual-trade-study method,
# re-exported under its original name so `weights_analysis.method` accepts either string. Not a
# copy: this always tracks `Physics_Based` exactly, since it imports the function directly rather
# than duplicating its body.
from ..Physics_Based.compute_operating_empty_weight import compute_operating_empty_weight
