
## @ingroup Library-Plots-Energy
# RCAIDE/Library/Plots/Energy/plot_l.py
# 
# 
# Created:  

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------  

from RCAIDE.Framework.Core import Units
from RCAIDE.Library.Plots.Common import set_axes, plot_style
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np 


# ----------------------------------------------------------------------------------------------------------------------
#  PLOTS
# ----------------------------------------------------------------------------------------------------------------------   
def plot_liquid_hydrogen_tank_properties(results,
                             save_figure = False,
                             show_legend = True,
                             save_filename = "Liquid_Hydrogen_Tank_Properties" ,
                             file_type = ".png", 
                             width = 11, height = 7):

    # get plotting style 
    ps = plot_style()  

    parameters = {
        'axes.labelsize': ps.axis_font_size,
        'xtick.labelsize': ps.axis_font_size,
        'ytick.labelsize': ps.axis_font_size,
        'axes.titlesize': ps.title_font_size
    }
    plt.rcParams.update(parameters)
     
    # get line colors for plots 
    line_colors = cm.inferno(np.linspace(0,0.9,len(results.segments)))      
         
    fig = plt.figure(save_filename)
    fig.set_size_inches(width, height) 

    axis_1 = plt.subplot(2,2,1)  # Mass
    axis_2 = plt.subplot(2,2,2)  # Temperatures
    axis_3 = plt.subplot(2,2,3)  # Volumes
    axis_4 = plt.subplot(2,2,4)  # Pressure
    
    for i, segment in enumerate(results.segments): 
        time = segment.conditions.frames.inertial.time[:, 0] / Units.min 

        for network in segment.analyses.energy.vehicle.networks: 
            for fuel_line in network.fuel_lines:
                for fuel_tank in fuel_line.fuel_tanks:

                    tank_conditions = segment.conditions.energy.fuel_lines[fuel_line.tag].fuel_tanks[fuel_tank.tag]

                    # Extract variables
                    tank_liquid_mass   = tank_conditions.mass[:, 0]
                    ullage_mass        = tank_conditions.ullage_mass[:, 0]
                    liquid_temp        = tank_conditions.liquid_temperature[:, 0]
                    ullage_temp        = tank_conditions.ullage_temperature[:, 0]
                    liquid_volume      = tank_conditions.liquid_volume[:, 0]/Units.gallons
                    ullage_volume      = tank_conditions.ullage_volume[:, 0]/Units.gallons
                    pressure           = tank_conditions.pressure[:, 0]
                    vent_rate          = tank_conditions.vent_rate[:, 0]
                    boil_off_rate = tank_conditions.boil_off_rate[:,0]

                    tank_label = f"{segment.tag}_{fuel_tank.tag}".replace("_"," ")

                    # --- Masses
                    axis_1.plot(time, tank_liquid_mass, color=line_colors[i],
                                linewidth=ps.line_width, label=f"{tank_label} liquid")
                    axis_1.plot(time, ullage_mass, color=line_colors[i],
                                linestyle="--", linewidth=ps.line_width, label=f"{tank_label} ullage")

                    # --- Temperatures
                    axis_2.plot(time, liquid_temp, color=line_colors[i],
                                linewidth=ps.line_width, label=f"{tank_label} liquid T")
                    axis_2.plot(time, ullage_temp, color=line_colors[i],
                                linestyle="--", linewidth=ps.line_width, label=f"{tank_label} ullage T")

                    # --- Volumes
                    axis_3.plot(time, liquid_volume, color=line_colors[i],
                                linewidth=ps.line_width, label=f"{tank_label} liquid V")
                    axis_3.plot(time, ullage_volume, color=line_colors[i],
                                linestyle="--", linewidth=ps.line_width, label=f"{tank_label} ullage V")

                    # --- Pressure
                    axis_4.plot(time, vent_rate, color=line_colors[i],
                                linewidth=ps.line_width, label=f"{tank_label} Vent Rate")
                    axis_4.plot(time, boil_off_rate, color='b',
                                linewidth=ps.line_width, label=f"{tank_label} Boiloff Rate")

    # Axis labels
    axis_1.set_ylabel("Mass (kg)")
    axis_2.set_ylabel("Temperature (K)")
    axis_3.set_ylabel("Volume (gal)")
    axis_4.set_ylabel("Mass Flow (kg/s)")

    for ax in [axis_1, axis_2, axis_3, axis_4]:
        ax.set_xlabel("Time (min)")
        set_axes(ax)

    # Legend
    if show_legend:
        leg = fig.legend(bbox_to_anchor=(0.5, 0.95), loc="upper center", ncol=4)

    # Adjust layout
    fig.tight_layout()
    fig.subplots_adjust(top=0.85)

    # Title
    fig.suptitle("Liquid Hydrogen Tank Properties")

    if save_figure:
        plt.savefig(save_filename + file_type)

    return fig