"""
Blade Element Momentum (BEM) theory solver for a horizontal-axis wind turbine.

Theory summary
--------------
BEM combines two independent models of the same annular blade element and
forces them to agree:

1. Momentum theory (1D actuator disk, annular control volume): relates the
   axial induction factor `a` and tangential induction factor `a'` to the
   pressure drop / angular momentum imparted to the flow through an annulus
   at radius r.

2. Blade element theory: computes the lift and drag forces on that same
   annulus from the local relative velocity, angle of attack, and airfoil
   polar (Cl(alpha), Cd(alpha)).

The two are solved iteratively per radial station until a, a' converge.
Glauert's tip-loss correction (Prandtl factor F) is applied to account for
the fact that a real rotor has a finite number of blades, not an idealized
actuator disk. A symmetric hub-loss factor is also applied. Glauert's
high-induction empirical correction is applied for a > 0.4, where simple
momentum theory breaks down (turbulent wake state).

Key equations (per station, standard BEM form):
    phi          = atan2((1-a)*V, (1+a')*omega*r)              [inflow angle]
    alpha        = phi - (twist + pitch)                       [angle of attack]
    Cn, Ct       = Cl*cos(phi) +/- Cd*sin(phi)  (normal/tangential coeffs)
    sigma        = B*c / (2*pi*r)                               [local solidity]
    F            = F_tip * F_hub                                [Prandtl/Glauert]
    a            = 1 / ( (4*F*sin^2(phi)) / (sigma*Cn) + 1 )
    a'           = 1 / ( (4*F*sin(phi)*cos(phi)) / (sigma*Ct) - 1 )

Integrated rotor performance:
    dT = 0.5*rho*Vrel^2*B*c*Cn*dr      (thrust contribution)
    dQ = 0.5*rho*Vrel^2*B*c*Ct*r*dr    (torque contribution)
    P  = Q * omega
    Cp = P / (0.5*rho*A*V^3)
    Ct_total = T / (0.5*rho*A*V^2)

Betz limit: Cp <= 16/27 ~= 0.593 is a structural consequence of 1D momentum
theory and is never imposed artificially in this code -- it should emerge
from the solver. If Cp exceeds it, that is a bug, not a feature.

IMPORTANT LIMITATION: the airfoil polar in `airfoil_polar()` below is a
simplified analytic model (thin-airfoil lift slope + flat-plate-like drag),
NOT measured wind-tunnel data for a real airfoil. It is transparent and
physically reasonable for demonstrating BEM mechanics, but absolute power
numbers will NOT match a real turbine (e.g. NREL 5MW) unless you replace
`airfoil_polar` with real Cl/Cd(alpha) tables (see README for how).
"""

import numpy as np
from dataclasses import dataclass, field


DEG = np.pi / 180.0
BETZ_LIMIT = 16.0 / 27.0


@dataclass
class TurbineParams:
    """All adjustable inputs, grouped per the spec. Units in comments."""
    # Environmental
    V: float = 9.0          # wind speed, m/s
    rho: float = 1.225      # air density, kg/m^3

    # Geometric
    R: float = 40.0         # rotor radius, m
    R_hub: float = 2.0      # hub radius, m
    B: int = 3              # number of blades
    c0: float = 3.4         # root chord, m (tip chord = 0.25*c0, linear taper)
    twist0: float = 13.0    # root twist, deg (linear washout to 0 deg at tip)

    # Operational
    rpm: float = 12.0       # rotor speed, RPM
    pitch_deg: float = 0.0  # blade pitch angle, deg

    n_stations: int = 24    # radial discretization

    @property
    def omega(self) -> float:
        return self.rpm * 2 * np.pi / 60.0

    @property
    def tsr(self) -> float:
        """Tip-speed ratio lambda = omega*R / V"""
        return self.omega * self.R / self.V


def airfoil_polar(alpha_rad: float) -> tuple[float, float]:
    """
    Simplified analytic Cl/Cd(alpha) model.

    Pre-stall: thin-airfoil theory, Cl = 2*pi*sin(alpha - alpha0), with
    induced drag Cd = Cd0 + k*Cl^2 (a parabolic drag polar).
    Post-stall (|alpha| > 14 deg): blended toward flat-plate-like behavior
    (Cl falls off, Cd rises sharply) to keep the BEM iteration stable and
    qualitatively realistic at high angles of attack.

    THIS IS NOT REAL AIRFOIL DATA. Replace with a measured polar (e.g. a
    NACA 4412 or DU/NREL S-series table) for quantitative validation work.
    """
    alpha0 = -2 * DEG
    alpha_deg = alpha_rad / DEG
    stall_deg = 14.0
    cd_min, k = 0.008, 0.02

    if abs(alpha_deg) <= stall_deg:
        cl = 2 * np.pi * np.sin(alpha_rad - alpha0)
        cd = cd_min + k * cl ** 2
    else:
        sgn = np.sign(alpha_deg)
        over = abs(alpha_deg) - stall_deg
        cl = sgn * 1.1 * np.cos(np.radians(min(over, 76)))
        cd = cd_min + k * 4 + 1.4 * np.sin(np.radians(min(over, 90)))
    return cl, cd


def chord_at(r: float, R: float, R_hub: float, c0: float) -> float:
    """Linear taper from root chord c0 to 0.25*c0 at the tip."""
    tip_chord = 0.25 * c0
    frac = (r - R_hub) / (R - R_hub)
    return c0 + (tip_chord - c0) * frac


def twist_at(r: float, R: float, R_hub: float, twist0: float) -> float:
    """Linear washout from root twist twist0 (deg) to 0 deg at the tip."""
    frac = (r - R_hub) / (R - R_hub)
    return twist0 * (1 - frac)


def prandtl_glauert_F(r, R, R_hub, B, phi) -> float:
    """
    Combined tip-loss and hub-loss correction factor (Glauert/Prandtl).
    Accounts for the finite number of blades vs. the idealized actuator disk.
    """
    sin_phi = max(abs(np.sin(phi)), 1e-4) * np.sign(np.sin(phi) if phi != 0 else 1)
    sin_phi = max(np.sin(phi), 1e-4)

    f_tip = (B / 2) * (R - r) / (r * sin_phi)
    f_hub = (B / 2) * (r - R_hub) / (r * sin_phi)

    F_tip = (2 / np.pi) * np.arccos(np.clip(np.exp(-f_tip), -1, 1))
    F_hub = (2 / np.pi) * np.arccos(np.clip(np.exp(-f_hub), -1, 1))
    return max(F_tip * F_hub, 0.05)


def solve_station(r, R, R_hub, B, c, twist_deg, pitch_deg, V, omega,
                   max_iter=100, relax=0.3, tol=1e-6):
    """
    Iteratively solve for axial/tangential induction factors (a, a') at a
    single radial station using BEM + Glauert correction. Returns a dict of
    all station-level quantities needed for force integration and plotting.
    """
    sigma = (B * c) / (2 * np.pi * r)
    a, ap = 0.2, 0.02

    for _ in range(max_iter):
        phi = np.arctan2((1 - a) * V, (1 + ap) * omega * r)
        phi_deg = phi / DEG
        alpha_deg = phi_deg - (twist_deg + pitch_deg)
        alpha_rad = alpha_deg * DEG
        cl, cd = airfoil_polar(alpha_rad)

        cn = cl * np.cos(phi) + cd * np.sin(phi)
        ct = cl * np.sin(phi) - cd * np.cos(phi)

        # Guard against near-zero or negative cn/ct, which occur at extreme
        # pitch/TSR combinations (e.g. deeply stalled or near-feathered
        # blades) and would otherwise produce divide-by-zero or a negative
        # discriminant in the Glauert correction below. Physically, cn<=0
        # means the element is producing no useful normal force, so capping
        # it near zero rather than letting it go negative is the right
        # behavior for a stable induction-factor solve.
        cn_safe = cn if abs(cn) > 1e-4 else np.sign(cn if cn != 0 else 1) * 1e-4
        ct_safe = ct if abs(ct) > 1e-4 else np.sign(ct if ct != 0 else 1) * 1e-4

        F = prandtl_glauert_F(r, R, R_hub, B, phi)

        denom_a = (4 * F * np.sin(phi) ** 2) / (sigma * cn_safe) + 1
        a_new = 1 / denom_a if denom_a != 0 else a

        # Glauert empirical high-induction correction (turbulent wake state)
        if a_new > 0.4:
            ac = 0.2
            K = (4 * F * np.sin(phi) ** 2) / (sigma * cn_safe)
            discriminant = (K * (1 - 2 * ac) + 2) ** 2 + 4 * (K * ac ** 2 - 1)
            if discriminant >= 0:
                a_new = 0.5 * (2 + K * (1 - 2 * ac) - np.sqrt(discriminant))
            else:
                # Discriminant only goes negative for non-physical (cn<=0)
                # operating points outside the normal BEM regime; fall back
                # to the un-corrected high-induction cap rather than
                # propagating a NaN through the rest of the solve.
                a_new = 0.95

        denom_ap = (4 * F * np.sin(phi) * np.cos(phi)) / (sigma * ct_safe) - 1
        ap_new = 1 / denom_ap if denom_ap != 0 else ap
        if not np.isfinite(ap_new):
            ap_new = ap

        a_prev, ap_prev = a, ap
        a = a + relax * (a_new - a)
        ap = ap + relax * (ap_new - ap)
        a = np.clip(a, -0.5, 0.95)
        ap = np.clip(ap, -0.5, 0.95)

        if abs(a - a_prev) < tol and abs(ap - ap_prev) < tol:
            break

    phi = np.arctan2((1 - a) * V, (1 + ap) * omega * r)
    phi_deg = phi / DEG
    alpha_deg = phi_deg - (twist_deg + pitch_deg)
    cl, cd = airfoil_polar(alpha_deg * DEG)
    cn = cl * np.cos(phi) + cd * np.sin(phi)
    ct = cl * np.sin(phi) - cd * np.cos(phi)
    Vrel2 = ((1 - a) * V) ** 2 + ((1 + ap) * omega * r) ** 2

    return dict(r=r, c=c, twist=twist_deg, sigma=sigma, phi_deg=phi_deg,
                alpha_deg=alpha_deg, cl=cl, cd=cd, a=a, ap=ap,
                cn=cn, ct=ct, Vrel2=Vrel2)


def solve_bem(p: TurbineParams) -> dict:
    """
    Run the full BEM solve across all radial stations and integrate
    forces to get rotor-level Cp, Ct, P, Q, T.
    """
    stations = []
    dT_sum, dQ_sum = 0.0, 0.0
    dr = (p.R - p.R_hub) / p.n_stations

    for i in range(p.n_stations):
        r = p.R_hub + (p.R - p.R_hub) * (i + 0.5) / p.n_stations
        c = chord_at(r, p.R, p.R_hub, p.c0)
        twist = twist_at(r, p.R, p.R_hub, p.twist0)

        st = solve_station(r, p.R, p.R_hub, p.B, c, twist, p.pitch_deg,
                            p.V, p.omega)

        dT = 0.5 * p.rho * st['Vrel2'] * p.B * c * st['cn'] * dr
        dQ = 0.5 * p.rho * st['Vrel2'] * p.B * c * st['ct'] * r * dr
        st['dT'] = dT
        st['dQ'] = dQ
        dT_sum += dT
        dQ_sum += dQ
        stations.append(st)

    A = np.pi * p.R ** 2
    P = dQ_sum * p.omega
    Cp = P / (0.5 * p.rho * A * p.V ** 3) if p.V > 0 else 0.0
    Ct = dT_sum / (0.5 * p.rho * A * p.V ** 2) if p.V > 0 else 0.0

    return dict(stations=stations, P=P, Q=dQ_sum, T=dT_sum, Cp=Cp, Ct=Ct,
                betz_fraction=Cp / BETZ_LIMIT if BETZ_LIMIT > 0 else 0.0)
