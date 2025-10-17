# RCAIDE/Methods/Powertrain/Sources/Fuel_Tanks/append_liquid_hydrogen_tank_unknown_and_residual.py
# 
# 
# Created:  Jul 2023, M. Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports 

# ----------------------------------------------------------------------------------------------------------------------
#  METHOD
# ----------------------------------------------------------------------------------------------------------------------    
def append_liquid_hydrogen_tank_unknown_and_residual(fuel_tank, segment):
    ones_row    = segment.state.ones_row

    segment.state.number_of_network_unknowns  += 6 
    segment.state.number_of_network_residuals += 6
    
    segment.state.unknowns.network[fuel_tank.tag + '_ullage_mass']                      = ones_row(1) * fuel_tank.ullage.mass_properties.mass 
    segment.state.unknowns.network[fuel_tank.tag + '_liquid_mass']                      = ones_row(1) * fuel_tank.fuel.mass_properties.mass 
    segment.state.unknowns.network[fuel_tank.tag + '_ullage_temperature']               = ones_row(1) * fuel_tank.ullage.temperature
    segment.state.unknowns.network[fuel_tank.tag + '_liquid_temperature']               = ones_row(1) * 20 # fuel_tank.fuel.temperature The defaulted value is not read here, I dont know why
    segment.state.unknowns.network[fuel_tank.tag + '_ullage_volume']                    = ones_row(1) * (fuel_tank.volume_properties.net_volume   -  fuel_tank.fuel.volume_properties.net_volume)
    segment.state.unknowns.network[fuel_tank.tag + '_liquid_volume']                    = ones_row(1) * fuel_tank.fuel.volume_properties.net_volume
               
           
    segment.state.residuals.network[fuel_tank.tag + '_ullage_mass']                     = ones_row(1) * 0
    segment.state.residuals.network[fuel_tank.tag + '_liquid_mass']                     = ones_row(1) * 0 
    segment.state.residuals.network[fuel_tank.tag + '_ullage_temperature']              = ones_row(1) * 0
    segment.state.residuals.network[fuel_tank.tag + '_liquid_temperature']              = ones_row(1) * 0 
    segment.state.residuals.network[fuel_tank.tag + '_ullage_volume']                   = ones_row(1) * 0
    segment.state.residuals.network[fuel_tank.tag + '_liquid_volume']                   = ones_row(1) * 0
 

    segment.state.unknowns_lower_bounds.network[fuel_tank.tag + '_ullage_mass']         = fuel_tank.ullage.mass_properties.mass * ones_row(1)
    segment.state.unknowns_lower_bounds.network[fuel_tank.tag + '_ullage_mass']         = fuel_tank.fuel.mass_properties.mass * ones_row(1) 
    segment.state.unknowns_lower_bounds.network[fuel_tank.tag + '_liquid_mass']         = 1e-6* ones_row(1)
    segment.state.unknowns_lower_bounds.network[fuel_tank.tag + '_liquid_mass']         = fuel_tank.fuel.mass_properties.mass* ones_row(1) 
    segment.state.unknowns_lower_bounds.network[fuel_tank.tag + '_ullage_temperature']  = 10 * ones_row(1)
    segment.state.unknowns_lower_bounds.network[fuel_tank.tag + '_ullage_temperature']  = 30 * ones_row(1)

    segment.state.unknowns_lower_bounds.network[fuel_tank.tag + '_liquid_temperature'] = 10 * ones_row(1) 
    segment.state.unknowns_lower_bounds.network[fuel_tank.tag + '_liquid_temperature'] = 30 * ones_row(1)  
    segment.state.unknowns_lower_bounds.network[fuel_tank.tag + '_ullage_volume']      = 1e-6* ones_row(1)
    segment.state.unknowns_lower_bounds.network[fuel_tank.tag + '_ullage_volume']      = fuel_tank.fuel.volume_properties.net_volume* ones_row(1) 
    segment.state.unknowns_lower_bounds.network[fuel_tank.tag + '_liquid_volume']      = 1e-6 * ones_row(1)
    segment.state.unknowns_lower_bounds.network[fuel_tank.tag + '_liquid_volume']      = fuel_tank.fuel.volume_properties.net_volume * ones_row(1)

    return
