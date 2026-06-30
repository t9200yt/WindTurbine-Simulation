# Wind turbine BEM simulator — formula reference

This document lists every physics equation used in the simulation, in the
order the solver applies them. It corresponds directly to the code in
`physics/bem.py` in the GitHub repository — each section below names the
function where that equation lives, so you can cross-reference.

---

## 1. Setup quantities

**Tip-speed ratio**

```
lambda = omega * R / V
```
- `omega` — rotor angular speed, rad/s (`omega = RPM * 2*pi / 60`)
- `R` — rotor radius, m
- `V` — wind speed, m/s

**Rotor swept area**

```
A = pi * R^2
```

---

## 2. Blade geometry (per radial station r)

The blade is discretized into N radial stations between the hub radius
`R_hub` and the tip `R`. At each station r:

**Chord (linear taper)**
```
c(r) = c0 + (0.25*c0 - c0) * (r - R_hub)/(R - R_hub)
```
`c0` is the root chord; the tip chord is fixed at 25% of the root chord.

**Twist (linear washout)**
```
beta(r) = beta0 * (1 - (r - R_hub)/(R - R_hub))
```
`beta0` is the root twist angle (deg); twist decreases linearly to 0° at
the tip.

**Local solidity**
```
sigma(r) = B * c(r) / (2 * pi * r)
```
`B` is the number of blades. Solidity is the fraction of the annular
swept area at radius r that is physically occupied by blade material.

---

## 3. Airfoil aerodynamics (lift and drag coefficients)

This is the **simplified analytic model** used in place of measured
wind-tunnel data (see the "Accuracy" note in the main README — this is
the single biggest source of deviation from a real turbine's numbers).

**Pre-stall (|alpha| <= 14°): thin-airfoil theory**
```
Cl(alpha) = 2*pi * sin(alpha - alpha0)
Cd(alpha) = Cd_min + k * Cl^2
```
- `alpha0 = -2°` (zero-lift angle of attack)
- `Cd_min = 0.008`, `k = 0.02` (parabolic induced-drag polar)

**Post-stall (|alpha| > 14°): blended toward flat-plate behavior**
```
Cl(alpha) = sign(alpha) * 1.1 * cos(min(|alpha|-14°, 76°))
Cd(alpha) = Cd_min + 4k + 1.4 * sin(min(|alpha|-14°, 90°))
```
This produces the characteristic post-stall behavior of falling lift and
sharply rising drag, without modeling a specific real airfoil's measured
stall characteristics.

---

## 4. Blade Element Momentum (BEM) iteration

This is the core of the simulation, solved independently at every radial
station, then integrated over the whole blade.

**Inflow angle** (angle between the rotor plane and the local relative wind)
```
phi = atan2( (1-a)*V , (1+a')*omega*r )
```

**Angle of attack**
```
alpha = phi - (beta(r) + theta)
```
`theta` is the blade pitch angle (deg), applied uniformly across the span.

**Force coefficients resolved into the rotor plane**
```
Cn = Cl*cos(phi) + Cd*sin(phi)     (normal to rotor plane -> thrust)
Ct = Cl*sin(phi) - Cd*cos(phi)     (tangential -> torque)
```

**Prandtl tip-loss factor**
```
f_tip = (B/2) * (R - r) / (r * sin(phi))
F_tip = (2/pi) * arccos( exp(-f_tip) )
```

**Prandtl hub-loss factor** (same form, mirrored at the hub)
```
f_hub = (B/2) * (r - R_hub) / (r * sin(phi))
F_hub = (2/pi) * arccos( exp(-f_hub) )
```

**Combined correction factor**
```
F = F_tip * F_hub
```
This accounts for the rotor having a finite number of discrete blades
rather than behaving like an idealized continuous actuator disk — induced
velocities fall off near the tip and hub, which F captures.

**Axial induction factor** (from momentum/blade-element balance)
```
a = 1 / ( 4*F*sin^2(phi) / (sigma*Cn) + 1 )
```

**Glauert high-induction correction** (applied when a > 0.4, where simple
momentum theory breaks down — the "turbulent wake state")
```
K = 4*F*sin^2(phi) / (sigma*Cn)
a_c = 0.2
a = 0.5 * ( 2 + K*(1-2*a_c) - sqrt( (K*(1-2*a_c)+2)^2 + 4*(K*a_c^2 - 1) ) )
```

**Tangential induction factor**
```
a' = 1 / ( 4*F*sin(phi)*cos(phi) / (sigma*Ct) - 1 )
```

These equations are solved **iteratively**: start with a guess for `a, a'`
(0.2 and 0.02), recompute `phi`, `alpha`, `Cl`, `Cd`, `Cn`, `Ct`, `F`, then
new `a, a'` values, repeat until they stop changing (convergence), typically
within 20–60 iterations per station.

---

## 5. Force and power integration

Once `a, a'` have converged at every station:

**Relative velocity squared at the station**
```
Vrel^2 = ((1-a)*V)^2 + ((1+a')*omega*r)^2
```

**Thrust contribution from this station** (annulus of width dr)
```
dT = 0.5 * rho * Vrel^2 * B * c(r) * Cn * dr
```

**Torque contribution from this station**
```
dQ = 0.5 * rho * Vrel^2 * B * c(r) * Ct * r * dr
```

Summed over all stations to get total thrust `T` and torque `Q`.

**Mechanical power**
```
P = Q * omega
```

**Power coefficient**
```
Cp = P / (0.5 * rho * A * V^3)
```

**Thrust coefficient**
```
Ct_total = T / (0.5 * rho * A * V^2)
```

---

## 6. The Betz limit

```
Cp_max = 16/27 ≈ 0.5926
```

This is the theoretical maximum fraction of kinetic energy in the wind
that *any* idealized rotor (actuator disk, no losses) can extract — it
falls directly out of 1D momentum theory and applies regardless of blade
design. It is not imposed as a rule in the code; it is a structural
property of the equations above. If a result ever exceeded it, that would
indicate a bug, not a better turbine.

Real turbines typically achieve Cp in the 0.35–0.48 range at their design
tip-speed ratio, well below Betz, due to tip losses, drag, and the
finite-blade effects the Glauert correction accounts for.

---

## 7. Optimizer objective

The optimizer (see `optimizer/optimize.py`) holds blade geometry and wind
speed fixed, and searches over pitch angle `theta` and tip-speed ratio
`lambda` to maximize:

```
maximize  P(theta, lambda)  subject to:
    theta  in [-10°, 20°]
    lambda in [2, 13]
```

solved via `scipy.optimize.differential_evolution` (a gradient-free global
optimizer), cross-checked against a brute-force grid search over the same
bounds. See the optimizer file's comments for why a gradient-based method
(e.g. `scipy.optimize.minimize`) was deliberately not used: the objective
embeds the iterative solve above, which is not smooth enough for reliable
gradient-based convergence.

---

## Where to go for real turbine data

Everything above is internally consistent and dimensionally correct BEM
theory. What it is *not* is validated against a specific real turbine,
because the airfoil polar (section 3) is analytic rather than measured.
To get numbers that match a real machine (e.g. the NREL 5MW reference
turbine), the equations in sections 1, 2, and 4–7 stay the same — only
section 3 (the airfoil polar) and the chord/twist functions in section 2
would need to be replaced with that turbine's actual published data.
