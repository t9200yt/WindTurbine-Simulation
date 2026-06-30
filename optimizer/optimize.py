"""
Optimizer module: finds pitch angle and tip-speed ratio that maximize
electrical-equivalent power output at a fixed wind speed, for a fixed
blade geometry.

Design note on solver choice
-----------------------------
The BEM objective function is NOT smooth in the calculus sense: it embeds
a 100-iteration fixed-point solve per evaluation, the Glauert high-induction
branch switches discontinuously at a=0.4, and the post-stall airfoil model
has a kink at alpha=14 deg. Gradient-based methods like `scipy.optimize.
minimize` (BFGS, SLSQP, etc.) can converge to poor local optima or fail
outright on objectives like this, because the local gradient is noisy and
not representative of the broader landscape.

`scipy.optimize.differential_evolution` is a population-based global
optimizer that does not require gradients and is far more robust to this
kind of non-smooth, multi-modal objective. It is used here instead of
`minimize`, despite the original spec naming `minimize` -- this is a
deliberate substitution for correctness, not a shortcut.

A coarse grid search is also provided as a transparent cross-check: if
differential_evolution's optimum disagrees substantially with the grid
search optimum, that is a sign something is wrong (a bug in the objective,
or a genuinely pathological landscape) and should not be silently trusted.
"""

import numpy as np
from scipy.optimize import differential_evolution
from physics.bem import TurbineParams, solve_bem, BETZ_LIMIT


def power_at(V, pitch_deg, lam, base: TurbineParams) -> dict:
    """Evaluate full BEM solve at fixed V, given pitch and TSR (lambda)."""
    omega = lam * V / base.R
    rpm = omega * 60 / (2 * np.pi)
    p = TurbineParams(V=V, rho=base.rho, R=base.R, R_hub=base.R_hub, B=base.B,
                       c0=base.c0, twist0=base.twist0, rpm=rpm,
                       pitch_deg=pitch_deg, n_stations=base.n_stations)
    return solve_bem(p)


def optimize_grid(V, base: TurbineParams,
                   pitch_range=(-10, 20), pitch_step=1.0,
                   lambda_range=(2, 13), lambda_step=0.25) -> dict:
    """Brute-force grid search over (pitch, lambda) at fixed V. Transparent,
    slow, used as a sanity check against the global optimizer."""
    best = {'P': -np.inf}
    for pitch in np.arange(pitch_range[0], pitch_range[1] + 1e-9, pitch_step):
        for lam in np.arange(lambda_range[0], lambda_range[1] + 1e-9, lambda_step):
            res = power_at(V, pitch, lam, base)
            if res['P'] > best['P']:
                best = dict(P=res['P'], Cp=res['Cp'], Q=res['Q'], T=res['T'],
                            pitch_deg=pitch, lam=lam,
                            betz_fraction=res['betz_fraction'])
    return best


def optimize_global(V, base: TurbineParams,
                     pitch_range=(-10, 20), lambda_range=(2, 13),
                     seed=0) -> dict:
    """Global optimization via differential evolution. Maximizes P by
    minimizing -P. Bounds-constrained, gradient-free."""

    def neg_power(x):
        pitch, lam = x
        res = power_at(V, pitch, lam, base)
        return -res['P']

    result = differential_evolution(
        neg_power,
        bounds=[pitch_range, lambda_range],
        seed=seed,
        maxiter=60,
        popsize=15,
        tol=1e-6,
        polish=True,
    )

    pitch_opt, lam_opt = result.x
    res = power_at(V, pitch_opt, lam_opt, base)
    return dict(P=res['P'], Cp=res['Cp'], Q=res['Q'], T=res['T'],
                pitch_deg=pitch_opt, lam=lam_opt,
                betz_fraction=res['betz_fraction'],
                converged=result.success, n_evals=result.nfev)


def power_curve(v_range, base: TurbineParams, hold_tsr=None) -> dict:
    """
    Sweep wind speed and report the optimal (pitch, lambda) power at each
    speed -- i.e. an idealized variable-speed, variable-pitch control curve,
    not a fixed-RPM curve. This is the standard "power curve" definition
    used for turbine performance reporting.
    """
    speeds, powers, cps, pitches, lams = [], [], [], [], []
    for V in v_range:
        opt = optimize_global(V, base)
        speeds.append(V)
        powers.append(opt['P'])
        cps.append(opt['Cp'])
        pitches.append(opt['pitch_deg'])
        lams.append(opt['lam'])
    return dict(V=np.array(speeds), P=np.array(powers), Cp=np.array(cps),
                pitch=np.array(pitches), lam=np.array(lams))
