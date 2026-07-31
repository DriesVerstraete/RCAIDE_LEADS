# RCAIDE/Library/Methods/Mass_Properties/Weight_Buildups/Electric/VTOL/Hydra/compute_rotor_weight.py
#
#

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
#  Constants
# ----------------------------------------------------------------------------------------------------------------------
# Direct port of Hydra's `src/Python/Stage_1/blade_wt_modelv2.py::blade_wt_modelv2`
# (`github.com/VahanaOpenSource/vtol_sizing`), confirmed live (`rotor_class.py:550`) — an
# independent physics/FEA-style blade structural sizing model, not a statistical regression like
# `NDARC`'s AFDD82/AFDD00 (`compute_rotor_weight.py` in the sibling folder). No tilt-factor term,
# unlike `NDARC`.
#
# Source is entirely SI internally (unlike NDARC's rotor/wing files, which are imperial-calibrated)
# — confirmed by every density/stress/modulus constant in the source (e.g. rho_al=2800 kg/m^3,
# sigma_skin=47e6 Pa). No unit conversion needed here.
#
# Dependency `spar_properties()` ported alongside (material lookup table, self-contained, no
# further dependencies). The source also references `hydra.crosssectiondata.init_materials()` — a
# module-level singleton init call — but this is dead weight for the actual computation: the
# airfoil section properties it would provide (`contour_pr`) are fully hardcoded inline in the
# source function itself (source's own comment: "values hardcoded because they are functions of
# airfoil section only"), so nothing in the real result depends on `crosssectiondata`. Not ported.

_CONTOUR_PR = {
    'Izzskin': 0.2952, 'Iyyskin': 0.00398, 'A_skin': 2.0392,
    'Atotal': 0.0822, 'mi1a_fill': 0.03456, 'mi1a_skin': 1.0052,
    'dsbyt': 2.0392, 'mi2a_fill': 0.006993,
}


def _spar_properties(material):
    """ Direct port of Hydra's `spar_properties.py::spar_properties` — material lookup table for
        rotor spar sizing. Source: http://asm.matweb.com/ (titanium/aluminum),
        http://www.acpsales.com/ (carbon variants), ZL values from the aerodesigntool repo. """
    spar = {}
    if material == 'titanium':
        spar = {'rho': 4500.0, 'sigma_y': 880.0e6, 'E': 114.0e9, 'G': 44.0e9, 'tau_y': 760.0e6}
    elif material == 'aluminum':
        spar = {'rho': 2700.0, 'sigma_y': 276.0e6, 'E': 71.0e9, 'G': 27.3e9, 'tau_y': 207.0e6}
    elif material == 'isolated_uniaxial_carbon':
        spar = {'rho': 1660.0, 'sigma_y': 450.0e6, 'E': 122.0e9, 'E2': 10.0e9, 'G': 5.0e9, 'tau_y': 70.0e9}
    elif material == 'isolated_090_carbon':
        spar = {'rho': 1600.0, 'sigma_y': 570.0e6, 'E': 70.0e9, 'E2': 70.0e9, 'G': 5.0e9, 'tau_y': 47.0e6}
    elif material == 'isolated_pm45_carbon':
        spar = {'rho': 1600.0, 'sigma_y': 110.0e6, 'E': 17.0e9, 'E2': 17.0e9, 'G': 47.0e9, 'tau_y': 210.0e6}
    elif material == 'uniaxial_carbon':
        spar = {'rho': 1660.0, 'sigma_y': 450.0e6, 'E': 122.0e9, 'G': 5.0e9}
    elif material == '090_carbon':
        spar = {'rho': 1660.0, 'sigma_y': 275.0e6, 'tau_y': 47.0e6, 'E': 70.0e9, 'G': 33.3e9}
    else:
        raise ValueError(f"_spar_properties: unknown spar material {material!r} "
                          "(expected 'titanium', 'aluminum', 'isolated_uniaxial_carbon', "
                          "'isolated_090_carbon', 'isolated_pm45_carbon', 'uniaxial_carbon', "
                          "or '090_carbon')")
    return spar


def compute_rotor_weight(radius, chord, omega, thrust, n_blade, n_rotor, material,
                          load_factor=1.0, precone_deg=3.0, rho_filler=52.0, tech_factor=0.5):
    """ Calculates rotor blade, hub, and pitch-actuator mass using Hydra's physics/FEA-style blade
        structural sizing model: an iterative spanwise centrifugal-force/bending-moment/spar-
        thickness sizing loop (5 stations), plus leading-edge protection, skin, filler, paint,
        glue, and root-fitting terms, capped off with a flat +20% for misc weights.

        Source:
            Hydra `blade_wt_modelv2.py::blade_wt_modelv2` + `spar_properties.py::spar_properties`.

        Inputs:
            radius         rotor radius                                              [m]
            chord          blade chord                                               [m]
            omega          rotor rotational speed                                    [rad/s]
            thrust         total rotor thrust for this rotor group                   [N]
            n_blade        number of blades per rotor                                [Unitless]
            n_rotor        number of rotors in this group                            [Unitless]
            material       spar material name (str, existing lookup-table behavior -- see
                            `_spar_properties`) OR a real `Solid` material instance (e.g.
                            `AS4_3502_Unidirectional_Carbon_Fiber`), 2026-07-30 -- only `.density`
                            (rho) and `.ultimate_tensile_strength` (sigma_y) are actually consumed
                            by this function's own sizing formula (confirmed by direct inspection
                            -- the `E`/`G`/`tau_y` entries in `_spar_properties`'s dict, and the
                            `web` variable below, are unused dead values in this specific port,
                            not just when a `Solid` instance is passed)                [str or Solid]
            load_factor    limit load factor (source applies a fixed 1.5 ultimate
                            factor on top of this internally)                        [Unitless]
            precone_deg    blade precone angle                                       [deg]
            rho_filler     honeycomb/foam filler density                             [kg/m^3]
            tech_factor    technology weight-scaling factor, BLADE ONLY -- hub/actuator
                            unaffected (default 0.5, decided 2026-07-30, see Finding T2,
                            00-decisions/2026-07-29-rotor-motor-weight-formula-comparison.md)   [Unitless]

        Outputs:
            weight:  dict with 'blades', 'hub', 'actuator', all summed over the whole
                     `n_rotor`-rotor group                                            [kg]
            total:   sum of 'blades'+'hub'+'actuator' — matches source's own second
                     return value, NOT divided by n_rotor (unlike NDARC's per-unit
                     convention)                                                      [kg]
    """
    if isinstance(material, str):
        spar = _spar_properties(material.lower())
    else:
        # real `Solid` material instance, 2026-07-30 -- see docstring: only rho/sigma_y consumed
        spar = {'rho': material.density, 'sigma_y': material.ultimate_tensile_strength}
    web = _spar_properties('090_carbon')  # unused below -- dead value, kept for source parity

    dts = 0.0005      # minimum thickness for manufacturing
    tbyc = 0.12       # blade thickness to chord ratio
    x_root = 0.1      # nondimensional span length of root fitting
    n_stations = 5    # number of spanwise stations

    betap = precone_deg * np.pi / 180.0   # precone angle, rad

    R = radius
    chord_ = chord
    nz = load_factor * 1.5
    Fz = thrust / n_blade

    r = np.linspace(x_root, 1.0, n_stations)
    b = 0.14 * chord_
    h = tbyc * chord_
    v_tip = omega * R

    c2 = chord_ * chord_
    c3 = chord_ * c2
    c4 = chord_ * c3

    # leading edge protection (Nickel plating)
    rho_lep = 8900.0
    l_lep = 0.34 * chord_
    a_lep = 2 * l_lep
    t_lep = 0.045 * 0.0254 * 0.5
    m_lep = a_lep * rho_lep * t_lep
    x_lep = 0.1 * chord_

    # skin sizing
    acs = _CONTOUR_PR['Atotal'] * c2
    tmin = 5e-4
    sigma_skin = 47.0e6
    rho_skin = 1660.0
    g_skin = 47.0e9

    cm = 0.02
    c1 = 1.2256 * v_tip * v_tip / 6.0 * c2 * cm * R
    ts_ref = max(nz * c1 / (2.0 * acs * sigma_skin), tmin)

    l_skin = _CONTOUR_PR['A_skin'] * chord_
    a_skin = l_skin * ts_ref
    m_skin = a_skin * rho_skin
    m_skin = m_skin + l_skin * 0.2   # lightning strike protection + bonding
    x_skin = 0.49 * chord_

    izz_skin = _CONTOUR_PR['Izzskin'] * c3 * ts_ref
    iyy_skin = _CONTOUR_PR['Iyyskin'] * c3 * ts_ref
    i_skin = (izz_skin + iyy_skin) * rho_skin

    # paint
    t_paint = 0.00015
    rho_paint = 1800.0
    t_ratio = t_paint / ts_ref
    rho_ratio = rho_paint / rho_skin
    m_paint = m_skin * t_ratio * rho_ratio
    x_paint = x_skin
    i_paint = i_skin * rho_ratio * t_ratio

    # glue
    t_glue = 2.54e-4
    rho_glue = 1800.0
    t_ratio_glue = t_glue / ts_ref
    rho_ratio_glue = rho_glue / rho_skin
    m_glue = t_glue * rho_glue * (l_skin + 0.2 * chord_)
    x_glue = x_skin
    i_glue = i_skin * rho_ratio_glue * t_ratio_glue

    # filler (foam/honeycomb)
    a_hc = chord_ * chord_ * tbyc * 0.4905
    m_hc = rho_filler * a_hc
    x_hc = 0.57 * chord_
    ixx_hc = _CONTOUR_PR['mi2a_fill'] * c4
    ixx_hc = ixx_hc - b * h * h * h / 12.0 - h * b * b * b / 12.0
    i_hc = ixx_hc * rho_filler

    # leading edge mass to move CG to quarter-chord
    m_total = m_lep + m_skin + m_hc + m_paint + m_glue
    x_cg = (m_lep * x_lep + m_skin * x_skin + m_hc * x_hc + m_paint * x_paint + m_glue * x_glue) / m_total

    x_le = 0.0
    m_le = m_total * (4.0 * x_cg / chord_ - 1.0)

    x_cg = (m_le * x_le + m_total * x_cg) / (m_le + m_total)
    m_total = m_total + m_le

    i_lep = m_lep * (0.25 * chord_ - x_lep) * (0.25 * chord_ - x_lep)
    i_le = m_le * 0.25 * 0.25 * c2
    i_total = i_hc + i_skin + i_lep + i_le + i_paint + i_glue

    # torsion natural frequency (computed but not enforced, matching source — informational only)
    ds = _CONTOUR_PR['dsbyt'] * chord_
    j_skin = 4 * acs * acs * ts_ref / ds
    gj_skin = g_skin * j_skin
    nutsq = 1.0 + np.pi * np.pi * 0.25 * gj_skin / (i_total * v_tip * v_tip)
    _nu_t = np.sqrt(nutsq)   # unused downstream, matches source (informational/diagnostic only)

    # spanwise sizing loop: centrifugal force, bending moment, spar thickness
    cf_out = 0.0
    a1 = m_total * v_tip * v_tip * 0.5
    a0 = spar['rho'] * (b + h) * v_tip * v_tip
    b1 = v_tip * v_tip * 0.5
    m_out = 0.0
    spar_mass = 0.0

    for i in range(n_stations - 1):
        xout = r[n_stations - i - 1]
        xin = r[n_stations - i - 2]
        dR = (xout - xin) * R
        c0 = b1 * (xout * xout - xin * xin)

        cfo = cf_out + a1 * (xout * xout - xin * xin)
        fbm = Fz * R * 0.25 * (xin**4 - 4 * xin + 3)

        d_rhs = (m_out - cf_out * dR * betap - m_total * dR * betap * c0 * 0.5) / (h * b)
        d_lhs = c0 * dR * betap * spar['rho'] * (b + h) / (b * h)
        rhs = cfo * 0.5 / (b + h) + abs(fbm) / (b * h) + d_rhs

        lhs = spar['sigma_y'] / nz - a0 * (xout * xout - xin * xin) * 0.5 / (b + h) + d_lhs
        t = max(rhs / lhs, tmin)

        vz = (1.0 - xin**3) * Fz
        tweb = vz / (h * sigma_skin / nz)

        a_spar = 2 * (b * t + h * tweb * 0.5)
        m = a_spar * spar['rho']
        spar_mass = spar_mass + m * dR
        m_out = m_out - cf_out * dR * betap - (m + m_total) * c0 * dR * betap * 0.5
        cf_out = cfo + m * c0

    m_total = m_total * R + spar_mass

    # root fitting
    sigma_al = 20e6
    rho_al = 2800.0
    r_root = 4 * h
    factor = 3.0

    eff_fbm = Fz * R * 0.25 * (r[0]**4 - 4 * r[0] + 3) + m_out
    rhs_root = cf_out / (2.0 * np.pi * r_root * sigma_al) + abs(eff_fbm) / (factor * np.pi * r_root * r_root * sigma_al)
    lhs_root = 1.0 - rho_al / sigma_al * omega * omega * (x_root * R) * 0.5
    t_root = rhs_root / lhs_root
    if t_root < 0.001:
        t_root = 0.001
    m_root = 2 * np.pi * r_root * t_root * (x_root * R) * rho_al
    m_total = m_total + m_root

    m_total = m_total * 1.2   # misc weight allowance

    mass_blade = m_total * n_blade * n_rotor

    mass_act = 1.47 * (omega / 228)**2 * (R / 0.75) * (chord_ / 0.082)**0.5
    mass_act = mass_act * n_rotor

    mass_hub = 4.84 * n_rotor * (R / 0.75)**0.5

    # 2026-07-30: tech_factor now scales BLADE ONLY, not hub/actuator. Real EASA TCDS anchors
    # (2026-07-29, 00-decisions/2026-07-29-rotor-motor-weight-formula-comparison.md, Finding T2)
    # showed the actuator term already tracks close to 1.0x real hardware, and imply hub is not
    # the driver of Hydra's heavy-outlier rotor mass either -- a blanket scale factor across all
    # three terms would incorrectly shrink hub/actuator along with blade. Default changed 1.0 ->
    # 0.5 at the same time (decided directly, not per-component-validated against real data --
    # rotor mass is a small fraction of total vehicle weight, so the residual imprecision this
    # leaves is accepted rather than pursuing a fully rigorous per-component fit).
    mass_blade = mass_blade * tech_factor
    total = mass_blade + mass_hub + mass_act

    weight = {'blades': mass_blade, 'hub': mass_hub, 'actuator': mass_act}

    return weight, total
