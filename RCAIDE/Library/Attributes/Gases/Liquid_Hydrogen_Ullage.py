# RCAIDE/Library/Attributes/Gases/Liquid_Hydrogen_Ullage.py
# 
# Created:  Mar 2024, S. Shekar

# ----------------------------------------------------------------------------------------------------------------------  
#  Imports
# ----------------------------------------------------------------------------------------------------------------------  
 
from .Gas import Gas  
from RCAIDE.Library.Components.Mass_Properties      import Mass_Properties

# ----------------------------------------------------------------------------------------------------------------------  
# Liquid_Hydrogen_Ullage Class
# ----------------------------------------------------------------------------------------------------------------------   

class Liquid_Hydrogen_Ullage(Gas):
    """
    A class representing carbon dioxide gas and its thermodynamic properties.

    Attributes
    ----------
    tag : str
        Identifier for the gas type ('Liquid_Hydrogen_Ullage')
  
    gas_specific_constant : float
        Specific gas constant in m²/s²-K
    composition : Container
        Chemical composition of the gas
            - CO2 : float
                Mass fraction of carbon dioxide (1.0 for pure CO2)

    Notes
    -----
    This class implements basic thermodynamic properties for carbon dioxide gas.
    All properties are for pure CO2 at standard conditions.
    
    **Definitions**
    
    'Specific Gas Constant'
        The individual gas constant for CO2, equal to the universal gas constant divided 
        by the molecular mass of CO2
    
    'Molecular Mass'
        The mass of one mole of CO2 molecules
    """
    def __defaults__(self):
        """This sets the default values.
        
            Assumptions:
                None
            
            Source:
                None
        """            
        self.tag                   ='Liquid_Hydrogen_Ullage'
        self.density               = 2.5 
        self.mass_properties       = Mass_Properties()