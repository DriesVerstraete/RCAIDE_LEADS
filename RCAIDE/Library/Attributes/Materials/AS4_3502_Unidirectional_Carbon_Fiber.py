# AS4_3502_Unidirectional_Carbon_Fiber.py
#
# Created: Jul 2026, D. Verstraete

#-------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------

from .Solid import Solid
from RCAIDE.Framework.Core import Units

#-------------------------------------------------------------------------------
# AS4/3502 Unidirectional Carbon Fiber Solid Class
#-------------------------------------------------------------------------------

class AS4_3502_Unidirectional_Carbon_Fiber(Solid):
    """
    A class representing a real, named, B-basis-tested unidirectional carbon fiber/epoxy tape
    system (AS4 12k fiber / 3502 epoxy resin), added 2026-07-30 as a rigorously-sourced
    alternative to `Unidirectional_Carbon_Fiber` (which is a generic, unnotched, room-temperature
    MatWeb median with no environmental or notch knockdown applied).

    Attributes
    ----------
    ultimate_tensile_strength : float
        0-degree (fiber-direction) compression B-basis allowable at 180F/wet, Pa (1000e6).
        Set to the COMPRESSION value, not tension, despite the field name (inherited from the
        `Solid` base class) -- a real bending/keel member is compression-critical on one flange,
        and 0-degree compression B-basis (1000 MPa) is lower than 0-degree tension B-basis
        (1379 MPa) for this system. Any caller using this field for a beam-bending sizing formula
        gets the physically controlling value; a caller that specifically wants tension should use
        `F1tu_180F_wet` documented below instead.
    ultimate_shear_strength : float
        In-plane shear B-basis allowable at 180F/wet, from a real +-45-degree laminate shear test
        (not a lamina-level estimate), Pa (81.4e6).
    ultimate_bearing_strength : float
        Not sourced this session -- no bearing table was located in the extracted chapter pages.
        Left at `Bidirectional_Carbon_Fiber`'s generic MatWeb value (600e6) as an interim
        placeholder. Needs a CMH-17 Volume 2/3 bearing-strength table (ASTM D5961) lookup,
        including its e/D and w/D geometry basis, before this is trustworthy.
    density : float
        Composite density, kg/m^3 (1570 -- nominal, real range 1560-1590).
    minimum_gage_thickness : float
        Nominal cured ply thickness, m (1.397e-4, i.e. 0.0055 in).
    youngs_modulus : float
        Not part of the `Solid` base interface -- added as a custom attribute for spar/beam
        natural-frequency sizing formulas. 0-deg tension modulus, mean (not B-basis), RT/dry,
        133e9 Pa (19.3 Msi).

    Notes
    -----
    Source: real B-basis test data for AS4 12k/3502 unidirectional tape, MIL-HDBK-17-2E-class
    "Fully Approved" statistical basis (Weibull/ANOVA-derived B-values, 30-40 specimens across
    5 batches per condition) -- extracted directly from
    `000-inbox/Carbon Fiber Composites_26_07_30_10_51_09.pdf`, Chapter 4, Section 4.2.8,
    Tables 4.2.8(a) [tension, 1-axis], 4.2.8(b) [tension, 1-axis, 250F/wet], 4.2.8(d)
    [compression, 1-axis], 4.2.8(g) [shear, 12-plane, +-45 laminate]. All strength values below
    are the "Normalized" column (specimen thickness + batch fiber volume normalized to 59% Vf).

    B-basis values by condition (ksi / MPa):
        Tension, 0-deg (F1tu):     RT/dry 205/1414,  180F/wet 200/1379,  250F/wet 191/1317
        Compression, 0-deg (F1cu): RT/dry 171/1179,  180F/wet 145/1000   (-65F/dry: Interim only, no B-value)
        Shear, +-45 laminate (F12su): RT/dry 13.4/92.4,  180F/wet 11.8/81.4,  250F/wet 10.3/71.0

    These are UNNOTCHED coupon/laminate B-basis allowables -- they do NOT include an open-hole/
    notch knockdown. Real open-hole (notched) B-basis data was searched for directly in both
    CMH-17-3H Volume 3 and its predecessor MIL-HDBK-17-3F (`000-inbox/MIL-HDBK-17-3F.pdf`,
    Chapter 4 "Building Block Approach") -- neither contains populated notched-strength VALUES for
    any material system; both only show "open hole tension/compression strength" as a row in a
    test-program-planning matrix (number of batches to run), never actual numbers. Real open-hole
    test data was only found in a different, non-handbook source: NASA/TM-2009-215900, "Hot/Wet
    Open Hole Compression Strength of Carbon/Epoxy Laminates" (real IM7/8552 test panels, layup
    `[+45,0,-45,0,90,0,0,90,0]S`, 5 of 9 plies at 0-deg -- directional, not pure UD, and not the
    same layup as this class's 100%-0-deg lamina data, so the ratio below is a rough
    order-of-magnitude proxy, not a same-basis notch-knockdown factor for AS4/3502 specifically):
        OHC, directional layup, RT dry, unconditioned:  73.7 ksi = 508 MPa
        OHC, directional layup, 220F/wet (0.7% moisture): 63.7 ksi = 439 MPa (-13.6% vs RT)
    Proxy retention ratio against this class's own 0-deg compression B-basis (180F/wet, 1000 MPa):
    439/1000 = ~44%. `Hydra`'s existing 450 MPa uniaxial bending constant sits almost exactly on
    this real open-hole hot/wet value (439-508 MPa) despite carrying no citation in RCAIDE --
    see `01-mission-profiles/00-decisions/2026-07-30-fuselage-weight-formula-comparison.md` for
    the full derivation and the Vahana MATLAB cross-check below.

    This value set is deliberately at 180F/wet (not RT/dry) to represent a real aircraft
    environmental design condition, consistent with the SF=1.5 conservatism already applied
    separately in RCAIDE's fuselage/wing weight-buildup formulas.

    Cross-check against the original Vahana study's own MATLAB source (`materials.m`,
    `999-software/vahanaTradeStudy/`, `github.com/VahanaOpenSource/vahanaTradeStudy`):
    `uni.stress = 450e6` Pa, explicitly commented "design ultimate tensile stress" -- i.e. the
    original Vahana team already intended this as a knocked-down design allowable, not a raw
    coupon value. RCAIDE's `Hydra` port carries this exact 450e6 constant forward unchanged;
    RCAIDE's own `Vahana`-named fuselage code (`Physics_Based/compute_fuselage_weight.py`)
    instead substitutes the generic, unnotched `Unidirectional_Carbon_Fiber` (1500 MPa) --
    i.e. RCAIDE's "Vahana" implementation has drifted from the real Vahana study's own material
    assumption, while "Hydra" preserved it faithfully.
    """

    def __defaults__(self):
        """Sets material properties at instantiation.

        Assumptions:
        None

        Source:
        MIL-HDBK-17-2E-class B-basis test data (see class docstring for exact tables/conditions).

        Inputs:
        N/A

        Outputs:
        N/A

        Properties Used:
        None
        """

        self.ultimate_tensile_strength  = 1000e6    * Units.Pa   # 0-deg compression B-basis, 180F/wet (controlling)
        self.ultimate_shear_strength    = 81.4e6    * Units.Pa   # +-45 laminate shear B-basis, 180F/wet
        self.ultimate_bearing_strength  = 600e6     * Units.Pa   # not sourced this session -- placeholder
        self.yield_tensile_strength     = 1000e6    * Units.Pa
        self.yield_shear_strength       = 81.4e6    * Units.Pa
        self.yield_bearing_strength     = 600e6     * Units.Pa
        self.minimum_gage_thickness     = 1.397e-4  * Units.m
        self.density                    = 1570.     * Units['kg/(m**3)']

        # Not part of the `Solid` base interface (no modulus field exists there) -- added here for
        # spar/beam natural-frequency sizing formulas (e.g. Hydra's wing/rotor spar models), which
        # need stiffness, not just strength. Real E1t mean (not B-basis -- modulus scatter is much
        # lower than strength scatter, mean is the normal reporting basis), RT/dry, Table 4.2.8(a).
        self.youngs_modulus              = 133e9     * Units.Pa   # 0-deg tension modulus, mean, RT/dry

        # Reference values not used by the Solid interface, kept for traceability
        self.F1tu_180F_wet              = 1379e6    * Units.Pa   # 0-deg tension B-basis, 180F/wet
        self.F1tu_RT_dry                = 1414e6    * Units.Pa   # 0-deg tension B-basis, RT/dry
        self.F1cu_RT_dry                = 1179e6    * Units.Pa   # 0-deg compression B-basis, RT/dry
        self.F12su_RT_dry               = 92.4e6    * Units.Pa   # +-45 shear B-basis, RT/dry

        # Open-hole (notched) reference data -- NOT from this same UD lamina/system, see docstring.
        # A different, directional (0-heavy) IM7/8552 laminate, real test data, NASA/TM-2009-215900.
        # Kept separate from the unnotched fields above -- apply explicitly, do not blend silently.
        self.OHC_directional_layup_RT_dry     = 508e6  * Units.Pa
        self.OHC_directional_layup_220F_wet   = 439e6  * Units.Pa
        self.OHC_proxy_knockdown_vs_F1cu_180F_wet = 0.44  # 439 MPa / 1000 MPa, order-of-magnitude only
