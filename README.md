# Wind Turbine BEM Simulator

An interactive simulator for a horizontal-axis wind turbine, built on 
Blade Element Momentum (BEM) theory — the same physics used in real 
turbine design.

Adjust wind speed, blade geometry (length, chord, twist, number of 
blades), and operating conditions (rotor speed, pitch angle) and see 
the effect on power output, torque, and efficiency in real time. 
Includes an animated rotor, performance curves, and an optimizer that 
finds the best pitch and rotor speed for a given wind speed.

## What it shows

- **Cp-λ curve** — how efficiently the rotor converts wind into power 
  across different tip-speed ratios, plotted against the theoretical 
  maximum (the Betz limit, ~59.3%)
- **Power curve** — power output across a range of wind speeds
- **Spanwise plots** — how chord, twist, angle of attack, and flow 
  induction vary along the blade, from hub to tip
- **Optimizer** — automatically searches for the pitch angle and 
  tip-speed ratio that maximize power at a chosen wind speed

## How to run it

```bash
pip install -r requirements.txt
streamlit run app.py
```

Then open the local link it prints in your browser.

## A note on accuracy

The aerodynamic model (lift/drag behavior of the blade) is a 
simplified analytical approximation, not measured wind-tunnel data 
from a real airfoil. It's accurate enough to demonstrate how the 
physics works and how each parameter affects performance, but the 
exact power numbers won't match a specific real-world turbine. See 
the comments in `physics/bem.py` for details.
