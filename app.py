"""
Streamlit dashboard for the HAWT BEM simulator.

Run with:
    streamlit run app.py

Layout: sidebar holds every adjustable input from spec section 3 (grouped
environmental / geometric / operational). Main area shows metrics, the
animated rotor, Cp-lambda curve, power curve, spanwise distributions, and
the optimizer panel.
"""

import streamlit as st
import numpy as np
from physics.bem import TurbineParams, solve_bem, BETZ_LIMIT
from viz.plots import rotor_animation_figure, cp_lambda_curve, power_curve_figure, spanwise_figures
from optimizer.optimize import optimize_global, optimize_grid, power_curve

st.set_page_config(page_title='HAWT BEM Simulator', layout='wide')
st.title('Wind turbine BEM simulator')
st.caption('Blade Element Momentum theory with Glauert tip/hub-loss correction. '
           'Airfoil polar is a simplified analytic model -- see README before treating '
           'absolute power numbers as quantitatively realistic for a specific turbine.')

with st.sidebar:
    st.header('Environmental')
    V = st.slider('Wind speed V (m/s)', 3.0, 25.0, 9.0, 0.5)
    rho = st.slider('Air density rho (kg/m^3)', 0.9, 1.4, 1.225, 0.005)

    st.header('Geometric')
    R = st.slider('Rotor radius R (m)', 10.0, 63.0, 40.0, 1.0)
    R_hub = st.slider('Hub radius R_hub (m)', 0.5, 5.0, 2.0, 0.1)
    B = st.slider('Number of blades B', 2, 4, 3, 1)
    c0 = st.slider('Root chord (m)', 1.0, 6.0, 3.4, 0.1)
    twist0 = st.slider('Root twist (deg)', 5.0, 25.0, 13.0, 0.5)

    st.header('Operational')
    rpm = st.slider('Rotor speed (RPM)', 2.0, 30.0, 12.0, 0.5)
    pitch_deg = st.slider('Pitch angle (deg)', -5.0, 25.0, 0.0, 0.5)

params = TurbineParams(V=V, rho=rho, R=R, R_hub=R_hub, B=int(B), c0=c0,
                        twist0=twist0, rpm=rpm, pitch_deg=pitch_deg)

st.sidebar.markdown(f'**Tip-speed ratio (lambda):** {params.tsr:.2f}')
st.sidebar.markdown(f'**Omega:** {params.omega:.2f} rad/s')

result = solve_bem(params)

col1, col2, col3, col4, col5 = st.columns(5)
col1.metric('Power', f"{result['P']/1000:.1f} kW")
col2.metric('Cp', f"{result['Cp']:.3f}")
col3.metric('% of Betz', f"{100*result['betz_fraction']:.1f}%")
col4.metric('Torque', f"{result['Q']:.0f} N·m")
col5.metric('Thrust', f"{result['T']:.0f} N")

tab_rotor, tab_cplambda, tab_power, tab_span, tab_opt = st.tabs(
    ['Rotor animation', 'Cp-lambda curve', 'Power curve', 'Spanwise distributions', 'Optimizer'])

with tab_rotor:
    st.plotly_chart(rotor_animation_figure(params), use_container_width=False)
    st.caption('Press Play. Rotation rate is illustrative, not real-time-accurate to RPM.')

with tab_cplambda:
    st.plotly_chart(cp_lambda_curve(params), use_container_width=True)
    st.caption(f'Betz limit = 16/27 = {BETZ_LIMIT:.4f}. The solver never exceeds this; '
               'if it did, that would indicate a bug.')

with tab_power:
    st.write('Idealized variable-speed/variable-pitch power curve '
             '(optimal pitch+lambda found at each wind speed). This sweep runs the optimizer '
             'at each point, so it is slower than the other tabs.')
    if st.button('Compute power curve (3-25 m/s)'):
        with st.spinner('Optimizing at each wind speed...'):
            pc = power_curve(np.arange(3, 26, 2), params)
        st.plotly_chart(power_curve_figure(pc['V'], pc['P'] / 1000), use_container_width=True)

with tab_span:
    figs = spanwise_figures(result['stations'])
    c1, c2 = st.columns(2)
    with c1:
        st.write('Chord c(r)')
        st.plotly_chart(figs['chord'], use_container_width=True)
        st.write('Angle of attack alpha(r)')
        st.plotly_chart(figs['aoa'], use_container_width=True)
    with c2:
        st.write('Twist beta(r)')
        st.plotly_chart(figs['twist'], use_container_width=True)
        st.write("Induction factors a, a'")
        st.plotly_chart(figs['induction'], use_container_width=True)

with tab_opt:
    st.write(f'Find pitch and tip-speed ratio that maximize power at the current wind speed '
              f'(**V = {V} m/s**), holding blade geometry fixed.')
    method = st.radio('Optimizer', ['Global (differential evolution)', 'Grid search (cross-check)'],
                       horizontal=True)
    if st.button('Run optimizer'):
        with st.spinner('Optimizing...'):
            if method.startswith('Global'):
                best = optimize_global(V, params)
            else:
                best = optimize_grid(V, params)
        st.success(
            f"**Optimal at V={V} m/s:** pitch = {best['pitch_deg']:.2f} deg, "
            f"lambda = {best['lam']:.2f}\n\n"
            f"**P\\* = {best['P']/1000:.2f} kW, Cp\\* = {best['Cp']:.4f} "
            f"({100*best['betz_fraction']:.1f}% of Betz limit)**"
        )
        st.caption('Run both methods and compare -- they should agree closely. '
                   'A large disagreement signals a problem with the objective landscape, not '
                   'something to paper over.')
