# blended_wing_body_aerodynamics_test.py

import RCAIDE
from RCAIDE.Framework.Core import Data, Units  
from RCAIDE.Library.Plots import *  
import numpy as  np 
import sys
import os
import matplotlib.pyplot as plt

sys.path.append(os.path.join( os.path.split(os.path.split(sys.path[0])[0])[0], 'Vehicles'))

# the analysis functions
from BWB    import vehicle_setup  ,  configs_setup

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main():
    
    vehicle  = vehicle_setup() 
    configs  = configs_setup(vehicle) 
    analyses = analyses_setup(configs)  
    mission  = mission_setup(analyses)
    missions = missions_setup(mission)  
    results  = missions.base_mission.evaluate() 

    vortex_distribution = results.segments.cruise.analyses.aerodynamics.settings.vortex_distribution
    plot_3d_vehicle_vlm_panelization(vortex_distribution=vortex_distribution,
                    save_filename               = "BWB_Top_View", 
                    show_wing_control_points    = False,  
                    show_figure                 =False)

    plot_3d_vehicle_vlm_panelization(vortex_distribution=vortex_distribution,
                    save_filename               = "BWB_Top_View",
                    show_wing_control_points    = True,  
                    show_figure                 =False)

    Cruise_CL        = results.segments.cruise.conditions.aerodynamics.coefficients.lift.total[2][0] 
    Cruise_CL_true   = 0.3841007724387933
    Cruise_CL_diff   = np.abs(Cruise_CL - Cruise_CL_true)
    
    
    plot_results(results)
    
    #print('Error: ',Cruise_CL_diff)
    #assert np.abs((Cruise_CL - Cruise_CL_true)/Cruise_CL_true) < 1e-6
    
    
    return 

# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def analyses_setup(configs):
    """Set up analyses for each of the different configurations."""

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # Build a base analysis for each configuration. Here the base analysis is always used, but
    # this can be modified if desired for other cases.
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle):
    """This is the baseline set of analyses to be used with this vehicle. Of these, the most
    commonly changed are the weights and aerodynamics methods."""

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Framework.Analyses.Vehicle()

    # ------------------------------------------------------------------
    #  Geometry
    # ------------------------------------------------------------------
    geometry = RCAIDE.Framework.Analyses.Geometry.Geometry()
    geometry.vehicle = vehicle
    geometry.settings.update_fuselage_properties = True
    geometry.settings.update_fuel_volume         = True
    geometry.settings.unique_geometry            = False
    analyses.append(geometry)
    

    # ------------------------------------------------------------------
    #  Weights
    weights = RCAIDE.Framework.Analyses.Weights.Conventional_BWB()
    weights.vehicle = vehicle 
    weights.settings.FLOPS.fidelity     = 'Complex'  
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method()
    aerodynamics.vehicle = vehicle   
    analyses.append(aerodynamics)
 
    # ------------------------------------------------------------------
    #  Energy
    energy = RCAIDE.Framework.Analyses.Energy.Energy()
    energy.vehicle = vehicle 
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = RCAIDE.Framework.Analyses.Planets.Earth()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    return analyses    
    
    

# ----------------------------------------------------------------------
#   Define the Mission
# ----------------------------------------------------------------------

def mission_setup(analyses):
    """This function defines the baseline mission that will be flown by the aircraft in order
    to compute performance."""

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'mission'
  
    Segments = RCAIDE.Framework.Mission.Segments 
    base_segment = Segments.Segment()

    # ------------------------------------------------------------------    
    #   Cruise Segment: Constant Speed Constant Altitude
    # ------------------------------------------------------------------    

    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise" 
    segment.analyses.extend( analyses.cruise ) 
    segment.altitude                                                 = 40000 * Units['ft']  
    segment.air_speed                                                = 450 * Units['knots']
    segment.distance                                                 = 2250 * Units.nmi 

    segment.state.numerics.mission_solver.verbose = True   
    
    # define flight dynamics to model             
    segment.flight_dynamics.force_x                                  = True  
    segment.flight_dynamics.force_z                                  = True     

    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['propulsor_1','propulsor_2', 'propulsor_3']] 
    segment.assigned_control_variables.body_angle.active             = True 

 
    segment.assigned_control_variables.hydrogen_tank_ullage_mass.active                          = True        
    segment.assigned_control_variables.hydrogen_tank_ullage_mass.assigned_hydrogen_tanks         = [['h2_aft_fuel_tank']]  
    segment.assigned_control_variables.hydrogen_tank_ullage_volume.active                        = True        
    segment.assigned_control_variables.hydrogen_tank_ullage_volume.assigned_hydrogen_tanks       = [['h2_aft_fuel_tank']]  
    segment.assigned_control_variables.hydrogen_tank_ullage_temperature.active                   = True        
    segment.assigned_control_variables.hydrogen_tank_ullage_temperature.assigned_hydrogen_tanks  = [['h2_aft_fuel_tank']]  
    segment.assigned_control_variables.hydrogen_tank_fuel_mass.active                            = True        
    segment.assigned_control_variables.hydrogen_tank_fuel_mass.assigned_hydrogen_tanks           = [['h2_aft_fuel_tank']]  
    segment.assigned_control_variables.hydrogen_tank_fuel_volume.active                          = True        
    segment.assigned_control_variables.hydrogen_tank_fuel_volume.assigned_hydrogen_tanks         = [['h2_aft_fuel_tank']]  
    segment.assigned_control_variables.hydrogen_tank_fuel_temperature.active                     = True        
    segment.assigned_control_variables.hydrogen_tank_fuel_temperature.assigned_hydrogen_tanks    = [['h2_aft_fuel_tank']]  

    mission.append_segment(segment) 

    return mission
 
def missions_setup(mission):
    """This allows multiple missions to be incorporated if desired, but only one is used here."""

    missions     = RCAIDE.Framework.Mission.Missions() 
    mission.tag  = 'base_mission'
    missions.append(mission)

    return missions

def plot_results(results):
    
    plot_liquid_hydrogen_tank_properties(results)
    return 
if __name__ == '__main__': 
    main()
    plt.show()