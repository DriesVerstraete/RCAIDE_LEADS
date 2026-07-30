# S2_449_SP381_Unidirectional_Glass_Fiber.py
#
# Created: Jul 2026, D. Verstraete

#-------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------

from .Solid import Solid
from RCAIDE.Framework.Core import Units

#-------------------------------------------------------------------------------
# S2-449/SP381 Unidirectional S-2 Glass Fiber Solid Class
#-------------------------------------------------------------------------------

class S2_449_SP381_Unidirectional_Glass_Fiber(Solid):
    """
    A class representing a real, named unidirectional S-2 glass fiber/epoxy tape system
    (Owens Corning S2-449 glass / 3M PR381 epoxy resin), added 2026-07-30 as a sourced glass-fiber
    alternative to the carbon classes, for cases where the higher strain-to-failure and damage
    tolerance of glass is preferred over carbon's higher stiffness.

    Attributes
    ----------
    ultimate_tensile_strength : float
        0-degree (fiber-direction) compression allowable, RT/dry, Pa (1158e6). Set to the
        COMPRESSION value despite the field name, for the same beam-bending-critical-flange
        reasoning as `AS4_3502_Unidirectional_Carbon_Fiber`. **Screening-class only (no B-basis
        available for this property at any condition) -- see Notes.**
    ultimate_shear_strength : float
        In-plane shear allowable, RT/dry, from a real +-45-degree laminate shear test, Pa
        (98.6e6). **Screening-class only.**
    ultimate_bearing_strength : float
        Not sourced this session. Left at `Bidirectional_Carbon_Fiber`'s generic MatWeb value
        (600e6) as an interim placeholder.
    density : float
        Composite density, kg/m^3 (1850 -- nominal, real range 1840-1970).
    minimum_gage_thickness : float
        Nominal cured ply thickness, m (2.286e-4, i.e. 0.009 in).

    Notes
    -----
    Source: `000-inbox/Glass Fiber Composites_26_07_30_10_51_33.pdf`, Chapter 6, Section 6.2.1,
    Tables 6.2.1(a) [tension, 1-axis], 6.2.1(d) [compression, 1-axis], 6.2.1(f) [shear, 12-plane,
    +-45 laminate], S2-449 43.5k/SP381 unidirectional tape.

    **Statistical basis is materially weaker than the carbon class.** Only 0-degree tension at
    RT/dry reached "Fully Approved" status (B-value = 198 ksi = 1365 MPa, 32 specimens/6 batches).
    Every other property/condition in this dataset -- including ALL of compression and ALL of
    shear at every temperature -- is "Screening" class only (2 batches, 9-12 specimens, no
    statistically-derived B-value, mean/CV only). The values used here (`ultimate_tensile_strength`,
    `ultimate_shear_strength`) are therefore Screening-class MEANS, not B-basis allowables -- a
    materially less conservative basis than the carbon class's B-values. Treat this class as a
    rough starting point, not a certifiable design allowable, until a fully-approved glass system
    is located (this specific dataset never reached full approval in this chapter) or the
    screening means are independently knocked down for statistical confidence.

    Values used (RT/dry, Screening class unless noted):
        Tension, 0-deg (F1tu):     B-value RT/dry 198 ksi = 1365 MPa (only B-basis in the dataset)
        Compression, 0-deg (F1cu): Screening mean RT/dry 168 ksi = 1158 MPa (no B-value at any condition)
        Shear, +-45 laminate (F12su): Screening mean RT/dry 14.3 ksi = 98.6 MPa,
                                        Screening mean 160F/wet 9.5 ksi = 65.5 MPa (no B-value at any condition)

    These are UNNOTCHED coupon/laminate values -- no open-hole/notch or environmental knockdown
    applied (RT/dry basis used throughout, since the hot/wet conditions in this dataset are
    Screening-class with very small sample sizes -- see
    `01-mission-profiles/00-decisions/2026-07-30-fuselage-weight-formula-comparison.md` for the
    carbon-fiber open-hole knockdown reference, NASA/TM-2009-215900).
    """

    def __defaults__(self):
        """Sets material properties at instantiation.

        Assumptions:
        None

        Source:
        MIL-HDBK-17-2E-class Screening/B-basis test data (see class docstring).

        Inputs:
        N/A

        Outputs:
        N/A

        Properties Used:
        None
        """

        self.ultimate_tensile_strength  = 1158e6    * Units.Pa   # 0-deg compression, Screening mean, RT/dry (controlling)
        self.ultimate_shear_strength    = 98.6e6    * Units.Pa   # +-45 laminate shear, Screening mean, RT/dry
        self.ultimate_bearing_strength  = 600e6     * Units.Pa   # not sourced this session -- placeholder
        self.yield_tensile_strength     = 1158e6    * Units.Pa
        self.yield_shear_strength       = 98.6e6    * Units.Pa
        self.yield_bearing_strength     = 600e6     * Units.Pa
        self.minimum_gage_thickness     = 2.286e-4  * Units.m
        self.density                    = 1850.     * Units['kg/(m**3)']

        # Reference values not used by the Solid interface, kept for traceability
        self.F1tu_RT_dry                = 1365e6    * Units.Pa   # 0-deg tension B-basis, RT/dry (the ONLY B-value here)
        self.F12su_160F_wet             = 65.5e6    * Units.Pa   # +-45 shear, Screening mean, 160F/wet
