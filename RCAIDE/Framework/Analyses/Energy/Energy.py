# RCAIDE/Framework/Analyses/Energy/Energy.py
# 
# 
# Created:  Jul 2023, M. Clarke 

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------  
# RCAIDE imports
from RCAIDE.Framework.Core     import Data
from RCAIDE.Framework.Analyses import Analysis 
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  ANALYSIS
# ----------------------------------------------------------------------------------------------------------------------  -Energy

# ----------------------------------------------------------------------------------------------------------------------
#  Energy Analysis
# ----------------------------------------------------------------------------------------------------------------------   
class Energy(Analysis):
    """ This is the base class for energy analyses.
    """
    def __defaults__(self):
        """This sets the default values and methods for the analysis.
            
        Assumptions:
            None

        Source:
            None 
        """        
        self.tag      = 'energy'
        self.vehicle  = Data()
        
    def evaluate(self,unknowns,state): 
        """Evaluate the thrust produced by the energy network.
    
        Assumptions:
            None

        Source:
            None

        Args:
            state (dict): flight conditions [-]

        Returns:
            results : results of the thrust evaluation method. 
        """ 
        # assumes only one network exists
        networks = self.vehicle.networks
        cg       = self.vehicle.mass_properties.center_of_gravity
        
        # Pack the unknowns to 
        if isinstance(unknowns,np.ndarray):
            state.unknowns.network.unpack_array(unknowns)

        networks.evaluate(state,cg)

        # Unpack Residuals
        residual_keys = list(state.residuals.network.keys())
        netowrk_res = Data()
        full_ures_vals = Data()
        for res in residual_keys:
            netowrk_res[res] = state.residuals.network[res]
            full_ures_vals[res] = netowrk_res[res] 

        return full_ures_vals.pack_array()
    