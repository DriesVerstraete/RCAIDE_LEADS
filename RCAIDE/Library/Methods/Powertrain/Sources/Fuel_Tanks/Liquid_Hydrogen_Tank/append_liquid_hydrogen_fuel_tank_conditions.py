# RCAIDE/Methods/Powertrain/Sources/Fuel_Tanks/append_liquid_hydrogen_fuel_tank_conditions.py
# 
# 
# Created:  Jul 2023, M. Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports
import  RCAIDE
from RCAIDE.Framework.Mission.Common     import   Conditions

# ----------------------------------------------------------------------------------------------------------------------
#  METHOD
# ----------------------------------------------------------------------------------------------------------------------  
def append_liquid_hydrogen_fuel_tank_conditions(tank, segment, distributor):
    """
    Appends initial conditions for liquid fuel tank component during later mission analysis.
    
    Parameters
    ----------
   
    which will be updated during mission analysis based on the fuel tank's performance
    and fuel consumption.
    
    See Also
    --------
    RCAIDE.Library.Components.Powertrain.Sources.Fuel_Tanks 
    """
    ones_row    = segment.state.ones_row
    
    distributor_conditions = segment.state.conditions.energy.fuel_lines[distributor.tag]
        
    distributor_conditions.fuel_tanks[tank.tag]                           = Conditions()  
    distributor_conditions.fuel_tanks[tank.tag].mass                      = 0 * ones_row(1)  
    distributor_conditions.fuel_tanks[tank.tag].mass_flow_rate            = 0 * ones_row(1)  
    distributor_conditions.fuel_tanks[tank.tag].surface_temperature       = 0 * ones_row(1)  
    distributor_conditions.fuel_tanks[tank.tag].boil_off_flow_rate        = 0 * ones_row(1)  
    distributor_conditions.fuel_tanks[tank.tag].ullage                    = 0 * ones_row(1)
    distributor_conditions.fuel_tanks[tank.tag].temperature               = 0 * ones_row(1)
    distributor_conditions.fuel_tanks[tank.tag].pressure                  = 0 * ones_row(1)
    distributor_conditions.fuel_tanks[tank.tag].volume_lh2                = 0 * ones_row(1)

    return 