# RCAIDE/Methods/Powertrain/Sources/Fuel_Tanks/append_liquid_hydrogen_fuel_tank_conditions.py
# 
# 
# Created:  Jul 2023, M. Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports
import  RCAIDE
from RCAIDE.Framework.Mission.Common     import   Conditions, Residuals, Unknowns

from copy import deepcopy

import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  METHOD
# ----------------------------------------------------------------------------------------------------------------------  
def append_liquid_hydrogen_tank_conditions(tank, segment, distributor):
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

    if tank.symmetric: # revisit later after things are working 
        tank.fuel.mass_properties.mass *= 0.5
        tank.ullage.mass_properties.mass *= 0.5
        tank.fuel.volume_properties.gross_volume *= 0.5
        tank.fuel.volume_properties.net_volume *= 0.5


    
    distributor_conditions = segment.state.conditions.energy.fuel_lines[distributor.tag]
        
    distributor_conditions.fuel_tanks[tank.tag]                           = Conditions()  
    distributor_conditions.fuel_tanks[tank.tag].mass_flow_rate            = ones_row(1) * 0
    distributor_conditions.fuel_tanks[tank.tag].boil_off_rate             = ones_row(1) * 0
    distributor_conditions.fuel_tanks[tank.tag].vent_rate                 = ones_row(1) * tank.vent_rate
    distributor_conditions.fuel_tanks[tank.tag].ullage_mass               = ones_row(1) * tank.ullage.mass_properties.mass
    distributor_conditions.fuel_tanks[tank.tag].mass                      = ones_row(1) * tank.fuel.mass_properties.mass
    distributor_conditions.fuel_tanks[tank.tag].ullage_temperature        = ones_row(1) * tank.ullage_temperature
    distributor_conditions.fuel_tanks[tank.tag].liquid_temperature        = ones_row(1) * tank.liquid_temperature
    distributor_conditions.fuel_tanks[tank.tag].ullage_volume             = ones_row(1) * (tank.fuel.volume_properties.gross_volume -  tank.fuel.volume_properties.net_volume)
    distributor_conditions.fuel_tanks[tank.tag].liquid_volume             = ones_row(1) * tank.fuel.volume_properties.net_volume
    distributor_conditions.fuel_tanks[tank.tag].pressure                  = ones_row(1) * 0
    distributor_conditions.fuel_tanks[tank.tag].secondary_fuel_flow_rate  = tank.secondary_fuel_flow_rate * ones_row(1) 
     
    if tank.symmetric:
        next_tag = tank.tag + "_symmetric"
        distributor_conditions.fuel_tanks[next_tag] = deepcopy(distributor_conditions.fuel_tanks[tank.tag])

    return 

def append_hydrogen_fuel_tank_segment_conditions(fuel_tank, segment, distributor): 

    # if type(distributor) == RCAIDE.Library.Components.Powertrain.Distributors.Electrical_Bus: 
    #     distributor_conditions = segment.state.conditions.energy.busses[distributor.tag]
    if  type(distributor) == RCAIDE.Library.Components.Powertrain.Distributors.Fuel_Line: 
        distributor_conditions = segment.state.conditions.energy.fuel_lines[distributor.tag]

    if segment.state.initials:  
        # if type(distributor) == RCAIDE.Library.Components.Powertrain.Distributors.Electrical_Bus: 
        #     distributor_initials = segment.state.initials.conditions.energy.busses[distributor.tag]
        if  type(distributor) == RCAIDE.Library.Components.Powertrain.Distributors.Fuel_Line: 
         distributor_initials = segment.state.initials.conditions.energy.fuel_lines[distributor.tag]
            
        distributor_conditions.fuel_tanks[fuel_tank.tag].ullage_mass[:,0]        = distributor_initials.fuel_tanks[fuel_tank.tag].ullage_mass[-1,0]
        distributor_conditions.fuel_tanks[fuel_tank.tag].mass[:,0]               = distributor_initials.fuel_tanks[fuel_tank.tag].mass[-1,0]
        distributor_conditions.fuel_tanks[fuel_tank.tag].ullage_temperature[:,0] = distributor_initials.fuel_tanks[fuel_tank.tag].ullage_temperature[-1,0]
        distributor_conditions.fuel_tanks[fuel_tank.tag].liquid_temperature[:,0] = distributor_initials.fuel_tanks[fuel_tank.tag].liquid_temperature[-1,0]
        distributor_conditions.fuel_tanks[fuel_tank.tag].ullage_volume[:,0]      = distributor_initials.fuel_tanks[fuel_tank.tag].ullage_volume[-1,0]
        distributor_conditions.fuel_tanks[fuel_tank.tag].liquid_volume[:,0]      = distributor_initials.fuel_tanks[fuel_tank.tag].liquid_volume[-1,0]
        distributor_conditions.fuel_tanks[fuel_tank.tag].pressure[:,0]           = distributor_initials.fuel_tanks[fuel_tank.tag].pressure[-1,0]

    return

def append_liquid_hydrogen_tank_residual_and_unknowns(fuel_tank, segment, distributor,network):
    ones_row    = segment.state.ones_row

    segment.state.number_of_network_unknowns  += 6 
    segment.state.number_of_network_residuals += 6

    distributor_unknowns   = segment.state.unknowns.network[network.tag].fuel_lines[distributor.tag]
    distributor_residuals  = segment.state.residuals.network[network.tag].fuel_lines[distributor.tag]
    distributor_lower_bounds   = segment.state.lower_bounds.network[network.tag].fuel_lines[distributor.tag]
    distributor_upper_bounds   = segment.state.upper_bounds.network[network.tag].fuel_lines[distributor.tag]

    distributor_unknowns[fuel_tank.tag] = Unknowns()
    distributor_residuals[fuel_tank.tag] = Residuals()
    
    distributor_unknowns[fuel_tank.tag].ullage_mass  = ones_row(1) * fuel_tank.ullage.mass_properties.mass 
    distributor_unknowns[fuel_tank.tag].mass  = ones_row(1) * fuel_tank.fuel.mass_properties.mass

    distributor_unknowns[fuel_tank.tag].ullage_temperature  = ones_row(1) * fuel_tank.ullage_temperature
    distributor_unknowns[fuel_tank.tag].liquid_temperature  = ones_row(1) * fuel_tank.liquid_temperature

    distributor_unknowns[fuel_tank.tag].ullage_volume  = ones_row(1) * (fuel_tank.fuel.volume_properties.gross_volume -  fuel_tank.fuel.volume_properties.net_volume)
    distributor_unknowns[fuel_tank.tag].liquid_volume  = ones_row(1) * fuel_tank.fuel.volume_properties.net_volume
    

    distributor_residuals[fuel_tank.tag].ullage_mass         = ones_row(1) * 0
    distributor_residuals[fuel_tank.tag].mass                = ones_row(1) * 0

    distributor_residuals[fuel_tank.tag].ullage_temperature  = ones_row(1) * 0
    distributor_residuals[fuel_tank.tag].liquid_temperature  = ones_row(1) * 0

    distributor_residuals[fuel_tank.tag].ullage_volume        = ones_row(1) * 0
    distributor_residuals[fuel_tank.tag].liquid_volume        = ones_row(1) * 0

    distributor_lower_bounds[fuel_tank.tag] = Conditions()
    distributor_upper_bounds[fuel_tank.tag] = Conditions()

    distributor_lower_bounds[fuel_tank.tag].ullage_mass = 0 * ones_row(1)
    distributor_upper_bounds[fuel_tank.tag].ullage_mass = np.inf * ones_row(1)#fuel_tank.ullage.mass_properties.mass * ones_row(1)
    
    distributor_lower_bounds[fuel_tank.tag].mass         = 0* ones_row(1)
    distributor_upper_bounds[fuel_tank.tag].mass         = np.inf * ones_row(1)

    distributor_lower_bounds[fuel_tank.tag].ullage_temperature = 5 * ones_row(1)
    distributor_upper_bounds[fuel_tank.tag].ullage_temperature = 35 * ones_row(1)

    distributor_lower_bounds[fuel_tank.tag].liquid_temperature = 5 * ones_row(1) 
    distributor_upper_bounds[fuel_tank.tag].liquid_temperature = 40 * ones_row(1) 

    distributor_lower_bounds[fuel_tank.tag].ullage_volume = 0* ones_row(1)
    distributor_upper_bounds[fuel_tank.tag].ullage_volume = np.inf * ones_row(1)#fuel_tank.fuel.volume_properties.net_volume* ones_row(1)

    distributor_lower_bounds[fuel_tank.tag].liquid_volume = 0 * ones_row(1)
    distributor_upper_bounds[fuel_tank.tag].liquid_volume = fuel_tank.fuel.volume_properties.net_volume * ones_row(1)

    return
