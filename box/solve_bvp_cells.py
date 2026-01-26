# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 20:32:00 2026

@author: simon
"""

import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

# %% Parameters
E = 1e9
nu = 0.3
G = E/(2*(1+nu))
g=9.8
A=0.0044
rho=1000
It = 10.513e-6
Iw = 1.64e-9
Ixx=1.14e-5
Iyy=4.18e-5
Ixy=1.41e-5
L = 2
M = 1000.0
plt_pause=0.5
# %% Differential equation for rotation
# GI_t*theta'' - EIw*theta'''' = m
def torsion_fun(x, y):
    # y[0]=theta, y[1]=theta', y[2]=theta'', y[3]=theta'''
    dydx = np.vstack((y[1], y[2], y[3], (G*It/E/Iw)*y[2]))
    return dydx

def torsion_bc(ya, yb):
    # ya = value at x=0, yb = values at x=L
    return np.array([
        ya[0],        # theta(0)=0
        ya[1],        # theta'(0)=0
        yb[2],        # theta''(L)=0
        G*It*yb[1] - E*Iw*yb[3] - M  # moment
    ])
# %% solve torsion
x = np.linspace(0, L, 100)
y_init = np.zeros((4, x.size))
sol = solve_bvp(torsion_fun, torsion_bc, x, y_init)
# %% plot rotation
fig1,ax=plt.subplots(num='rotation',clear=True)
ax.plot(sol.x, 180/np.pi*sol.y[0], label="rotation θ(x)")
ax.set_xlabel("x [m]")
ax.set_ylabel("θ(x) [deg]")
ax.legend()
plt.pause(plt_pause)
# %% plot moments
T = G*It*sol.y[1] - E*Iw*sol.y[3]
fig2,ax=plt.subplots(num='moment',clear=True)
ax.plot(sol.x, T, label="Total moment T(x)")
ax.set_xlabel("x [m]")
ax.set_ylabel("Moment [Nm]")
T = G*It*sol.y[1]
ax.plot(sol.x, T, label="Moment from It")
T = - E*Iw*sol.y[3]
ax.plot(sol.x, T, label="Moment from Iw")
ax.legend()
plt.pause(plt_pause)
# %% beam bending
def beam_bending_with_elastic_supports(L, EI, bc,
                               point_loads=None, q_func=None,
                               n_points=50):
    if q_func is None:
        q_func = lambda x: np.zeros_like(x)

    def q_total(x):
        q = q_func(x)
        if point_loads:
            for xp, P in point_loads:
                q += P * np.exp(-((x-xp)**2)/(2*1e-4)) / np.sqrt(2*np.pi*1e-4)
        return q

    def bending_fun(x, y):
        return np.vstack((y[1], y[2], y[3], q_total(x)/EI))
    
    xs = np.linspace(0, L, n_points)
    ya = np.zeros((4,xs.size))
    sol=solve_bvp(bending_fun, bc, xs, ya)
    reactions = []
    total_load = np.trapezoid(q_total(xs), xs)
    if point_loads:
        total_load += sum(P for _, P in point_loads)
    total_reaction = sum(R for _, R in reactions)
    equilibrium_error = total_load + total_reaction

    return (sol,reactions, xs, q_total,
            total_load, equilibrium_error)

def run_example(bc,point_loads,distributed_load):
    (res, reactions, xs, q_total, total_load,
     equilibrium_error) = (
        beam_bending_with_elastic_supports(
        L, E*Ixx, bc,
        point_loads=point_loads,
        q_func=distributed_load
    ))
    if not res.success:
        raise Exception(f"solution failed, {res.message}")
    y_plot, M_plot, V_plot = [], [], []
    for x in xs:
        y, dy, M, V = res.sol(x)
        y_plot.append(y)
        M_plot.append(M)
        V_plot.append(V)
    
    # Deflection
    fig1,ax1=plt.subplots(num='Deflection',clear=True)
    ax1.plot(xs, y_plot, label="Deflection")
    if supports:
        for xm, _ in supports:
            ax1.axvline(x=xm, color="blue", linestyle="--", 
                        alpha=0.5, label="Support" 
                        if xm==supports[0][0] else "")
    if point_loads:
        for xp, P in point_loads:
            ax1.axvline(x=xp, color="red", linestyle=":", 
                        alpha=0.7, label="Point load" 
                        if xp==point_loads[0][0] else "")
    ax1.set_ylabel("y(x)")
    ax1.legend()
    plt.pause(plt_pause)
    
    # Annotate reactions
    for xm, R in reactions:
        plt.annotate(f"R={R:.1f}", xy=(xm,0), xytext=(xm,0.05),
                     arrowprops=dict(facecolor='black', shrink=0.05),
                     ha='center')
    
    # Moment
    fig2,ax2=plt.subplots(num='Moment',clear=True)
    ax2.plot(xs, M_plot, label="Moment")
    if supports:
        for xm, _ in supports:
            ax2.axvline(x=xm, color="blue", linestyle="--", alpha=0.5)
    if point_loads:
        for xp, P in point_loads:
            ax2.axvline(x=xp, color="red", linestyle=":", alpha=0.7)
    ax2.set_ylabel("M(x)")
    ax2.set_xlabel("x")
    plt.pause(plt_pause)
    
    # Shear
    fig3,ax3=plt.subplots(num='Shear',clear=True)
    ax3.plot(xs, V_plot, label="Shear")
    if supports:
        for xm, _ in supports:
            ax3.axvline(x=xm, color="blue", linestyle="--", alpha=0.5)
    if point_loads:
        for xp, P in point_loads:
            ax3.axvline(x=xp, color="red", linestyle=":", alpha=0.7)
    ax3.set_ylabel("V(x)")
    ax3.set_xlabel("x")
    plt.pause(plt_pause)
    
    fig4,ax4=plt.subplots(num='Load',clear=True)
    ax4.plot(xs, q_total(xs), color="green", label="Load q(x)")
    if point_loads:
        for xp, P in point_loads:
            ax4.axvline(x=xp, color="red", linestyle=":", alpha=0.7)
    ax4.set_ylabel("q(x)")
    ax4.set_xlabel("x")
    ax4.legend()
    plt.pause(plt_pause)
    return res
# %% Clamped beam with own weight
def bc(ya,yb):
    return np.array([
        ya[0],# displacement(0)=0
        ya[1],# rotation'(0)=0
        yb[2],# moment''(L)=0
        yb[3] # shear'''(L)=0
    ])
point_loads=None
distributed_load = lambda x: -A*rho*g*np.ones_like(x)
res1=run_example(bc,point_loads,distributed_load)
# %% Clamped beam with load at free end
def bc(ya,yb):
    return np.array([
        ya[0],# displacement(0)=0
        ya[1],# rotation'(0)=0
        yb[2],# moment''(L)=0
        yb[3] # shear'''(L)=0
    ])
point_loads = [(1*L, -0.5*L*A*rho*g)]
distributed_load = None
res2=run_example(bc,point_loads,distributed_load)
# %% No loads and no supports
def bc(ya,yb):
    return np.array([
        ya[0],# displacement(0)=0
        ya[1],# rotation'(0)=0
        yb[0],# displacement(L)=0
        yb[1] # rotation'(L)=0
    ])
supports = None
point_loads = None
distributed_load = None
res3=run_example(bc,point_loads,distributed_load)
# %% Own weight with simple supports
def bc(ya,yb):
    return np.array([
        ya[0],# displacement(0)=0
        ya[2],# moment''(0)=0
        yb[0],# displacement''(L)=0
        yb[2] # moment''(L)=0
    ])
supports = None
point_loads = None
distributed_load = lambda x: -A*rho*g*np.ones_like(x)
run_example(bc,point_loads,distributed_load)

