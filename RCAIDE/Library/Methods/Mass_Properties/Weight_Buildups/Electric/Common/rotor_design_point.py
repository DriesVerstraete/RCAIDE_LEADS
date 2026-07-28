# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/Common/rotor_design_point.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
#  Rotor design-point helpers
# ----------------------------------------------------------------------------------------------------------------------
def max_design_field(rotor, field):
    """Returns the more severe of rotor.hover.<field> and rotor.cruise.<field>, treating the
    structural design load as the worst case across whichever design points the rotor actually
    has. A plain Propeller only ever has .cruise; Lift_Rotor/Prop_Rotor/base Rotor have both. Falls
    back to whichever single value is set/non-None if only one design point exists.
    """
    values = []
    hover = getattr(rotor, 'hover', None)
    if hover is not None:
        value = getattr(hover, field, None)
        if value is not None:
            values.append(value)
    cruise = getattr(rotor, 'cruise', None)
    if cruise is not None:
        value = getattr(cruise, field, None)
        if value is not None:
            values.append(value)
    if not values:
        raise ValueError(
            f"max_design_field: rotor '{rotor.tag}' has neither hover.{field} nor cruise.{field} set."
        )
    return max(values)
