"""
Visualization module. Pure functions: take BEM results in, return Plotly
figures out. No Streamlit calls here -- keeps this testable and reusable
outside the dashboard (e.g. for the README's static validation plots).
"""

import numpy as np
import plotly.graph_objects as go
from physics.bem import TurbineParams, solve_bem, BETZ_LIMIT


def rotor_animation_figure(p: TurbineParams, n_frames=24):
    """3-blade (or B-blade) rotor animation, front-on view, rotating at the
    angular rate implied by the current RPM (visual rate, not real-time)."""
    R = p.R
    frames = []
    theta0 = np.linspace(0, 2 * np.pi, n_frames, endpoint=False)

    def blade_traces(offset):
        traces = []
        for i in range(p.B):
            ang = offset + i * (2 * np.pi / p.B)
            tip_x, tip_y = R * np.sin(ang), R * np.cos(ang)
            w = R * 0.05
            w1x, w1y = w * np.sin(ang + np.pi / 2), w * np.cos(ang + np.pi / 2)
            w2x, w2y = w * np.sin(ang - np.pi / 2), w * np.cos(ang - np.pi / 2)
            traces.append(go.Scatter(
                x=[w1x, tip_x, w2x, w1x], y=[w1y, tip_y, w2y, w1y],
                fill='toself', fillcolor='#1baf7a', line=dict(color='#0f6e56', width=1),
                showlegend=False, hoverinfo='skip'
            ))
        return traces

    fig = go.Figure(
        data=blade_traces(0),
        layout=go.Layout(
            xaxis=dict(range=[-R * 1.15, R * 1.15], visible=False, scaleanchor='y'),
            yaxis=dict(range=[-R * 1.15, R * 1.15], visible=False),
            width=420, height=420, margin=dict(l=10, r=10, t=10, b=10),
            updatemenus=[dict(type='buttons', showactive=False,
                               buttons=[dict(label='Play', method='animate',
                                             args=[None, dict(frame=dict(duration=60, redraw=True),
                                                               fromcurrent=True)])])]
        ),
        frames=[go.Frame(data=blade_traces(t)) for t in theta0]
    )
    fig.add_shape(type='circle', x0=-R * 1.05, y0=-R * 1.05, x1=R * 1.05, y1=R * 1.05,
                   line=dict(color='#d3d1c7', dash='dot'))
    return fig


def cp_lambda_curve(base: TurbineParams, lam_range=np.arange(2, 14, 0.5)):
    cps = []
    for lam in lam_range:
        omega = lam * base.V / base.R
        rpm = omega * 60 / (2 * np.pi)
        p = TurbineParams(**{**base.__dict__, 'rpm': rpm})
        res = solve_bem(p)
        cps.append(min(res['Cp'], BETZ_LIMIT))

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=lam_range, y=cps, mode='lines', name='Cp(lambda)',
                              line=dict(color='#2a78d6', width=2)))
    fig.add_hline(y=BETZ_LIMIT, line=dict(color='#e34948', dash='dash'),
                  annotation_text='Betz limit (16/27)')
    fig.update_layout(xaxis_title='Tip-speed ratio (lambda)', yaxis_title='Cp',
                       height=320, margin=dict(l=40, r=20, t=20, b=40))
    return fig


def power_curve_figure(speeds, powers_kw):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=speeds, y=powers_kw, mode='lines', fill='tozeroy',
                              line=dict(color='#1baf7a', width=2)))
    fig.update_layout(xaxis_title='Wind speed V (m/s)', yaxis_title='Power (kW)',
                       height=320, margin=dict(l=40, r=20, t=20, b=40))
    return fig


def spanwise_figures(stations: list[dict]):
    r = [s['r'] for s in stations]
    figs = {}

    def line(y, name, color):
        f = go.Figure()
        f.add_trace(go.Scatter(x=r, y=y, mode='lines', name=name,
                                line=dict(color=color, width=2)))
        f.update_layout(height=220, margin=dict(l=40, r=10, t=10, b=30),
                         xaxis_title='r (m)')
        return f

    figs['chord'] = line([s['c'] for s in stations], 'chord', '#2a78d6')
    figs['twist'] = line([s['twist'] for s in stations], 'twist', '#1baf7a')
    figs['aoa'] = line([s['alpha_deg'] for s in stations], 'alpha', '#eda100')

    f = go.Figure()
    f.add_trace(go.Scatter(x=r, y=[s['a'] for s in stations], name='a',
                            line=dict(color='#2a78d6', width=2)))
    f.add_trace(go.Scatter(x=r, y=[s['ap'] for s in stations], name="a'",
                            line=dict(color='#e34948', width=2)))
    f.update_layout(height=220, margin=dict(l=40, r=10, t=10, b=30), xaxis_title='r (m)')
    figs['induction'] = f

    return figs
