# RCAIDE/Methods/Powertrain/Converters/Fuel_Cells/Larminie_Model/append_fuel_cell_unknown_and_residual.py
# 
# 
# Created:  Oct 2025, M. Clarke 

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports
import RCAIDE
# ----------------------------------------------------------------------------------------------------------------------
#  append_fuel_cell_unknown_and_residual
# ----------------------------------------------------------------------------------------------------------------------
def append_fuel_cell_unknown_and_residual(fuel_cell,segment): 

    # compute ambient conditions
    atmosphere    = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    alt           = -segment.conditions.frames.inertial.position_vector[:,2] 
    if segment.temperature_deviation != None:
        temp_dev = segment.temperature_deviation    
    atmo_data    = atmosphere.compute_values(altitude = alt,temperature_deviation=temp_dev)  
    ones_row    = segment.state.ones_row

    #if segment.initial_battery_conditions.cell_temperature is not None:
        #stack_temperature  = segment.battery_cell_temperature  
    #else:
        #stack_temperature = atmo_data.temperature[0,0]
        
    # Conditions for recharging fuel_cell   
    segment.state.unknowns.network[ fuel_cell.tag + '_current_density']                =  1 * ones_row(1)  
    #segment.state.unknowns.network[ fuel_cell.tag + '_stack_temperature']              =  stack_temperature * ones_row(1)  
    segment.state.residuals.network[ fuel_cell.tag + '_power']                         =  0* ones_row(1) 
    #segment.state.residuals.network[ fuel_cell.tag + '_stack_temperature']             =  0* ones_row(1) 
    segment.state.unknowns_upper_bounds.network[ fuel_cell.tag + '_current_density']   =  2 * ones_row(1) 
    segment.state.unknowns_lower_bounds.network[ fuel_cell.tag + '_current_density']   =  0 * ones_row(1) 
    #segment.state.unknowns_upper_bounds.network[ fuel_cell.tag + '_stack_temperature'] =  350 * ones_row(1) 
    #segment.state.unknowns_lower_bounds.network[ fuel_cell.tag + '_stack_temperature'] =  200 * ones_row(1) 
    segment.state.number_of_network_unknowns  += 1
    segment.state.number_of_network_residuals += 1          
    
    return 
