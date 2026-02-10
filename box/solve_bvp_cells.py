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
res={}
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
# %% bc_multi sample
"""
By AI
Got it ✅ — here’s a general-purpose multi-region BVP solver template in
Python using scipy.integrate.solve_bvp that works for any number of regions,
each with its own ODE definition and interface continuity conditions.
This template is fully runnable, handles arbitrary region boundaries,
and enforces continuity of yyy and y′y'y′ at each interface.
# ============================================================
# 1. Define your problem
# ============================================================

# Example: y'' = -k^2 * y in each region, with different k
region_boundaries = [0, 1, 2, 3]  # 3 regions: [0,1], [1,2], [2,3]
region_k = [1.0, 2.0, 0.5]        # k for each region

# ODE system: y[0] = y, y[1] = y'
def fun(x, y):
    dydx = np.zeros_like(y)
    dydx[0] = y[1]
    # Determine which region each x belongs to
    k_values = np.zeros_like(x)
    for i in range(len(region_k)):
        mask = (x >= region_boundaries[i]) & (x <= region_boundaries[i+1])
        k_values[mask] = region_k[i]
    dydx[1] = -(k_values**2) * y[0]
    return dydx

# ============================================================
# 2. Boundary conditions
# ============================================================

# General multipoint BC function
def bc_multi(*Y):
    # Y[0] = values at start of region 1
    # Y[1] = values at end of region 1 (= start of region 2)
    # ...
    #Y[-1] = values at end of last region
    conditions = []
    # Left boundary condition: y(0) = 0
    conditions.append(Y[0][0] - 0.0)
    # Right boundary condition: y(end) = 0
    conditions.append(Y[-1][0] - 0.0)
    # Continuity conditions at internal boundaries
    for i in range(1, len(Y)-1):
        # Continuity of y
        conditions.append(Y[i][0] - Y[i+1][0])
        # Continuity of y'
        conditions.append(Y[i][1] - Y[i+1][1])
    return np.array(conditions)
# ============================================================
# 3. Build mesh and initial guess
# ============================================================

# Create mesh with duplicated interface points
mesh_parts = []
for i in range(len(region_boundaries)-1):
    part = np.linspace(region_boundaries[i], region_boundaries[i+1], 5)
    if i > 0:
        part = part[1:]  # avoid duplicate except first
    mesh_parts.append(part)
x_mesh = np.hstack(mesh_parts)

# Initial guess: y=0, y'=0
y_guess = np.zeros((2, x_mesh.size))

# ============================================================
# 4. Wrap BC function for solve_bvp
# ============================================================
def bc_wrapper(Ya, Yb):
    #    Ya: shape (2, n_regions) — values at left ends
    #    Yb: shape (2, n_regions) — values at right ends
    # Combine into sequence: start, interface1, interface2, ..., end
    Y_points = [Ya[:, 0]]
    for i in range(Ya.shape[1]):
        Y_points.append(Yb[:, i])
    return bc_multi(*Y_points)
# ============================================================
# 5. Solve
# ============================================================

sol = solve_bvp(fun, bc_wrapper, x_mesh, y_guess)

if sol.status != 0:
    print("Warning: Solver did not converge.")

# ============================================================
# 6. Plot
# ============================================================

x_plot = np.linspace(region_boundaries[0], region_boundaries[-1], 300)
y_plot = sol.sol(x_plot)[0]

plt.figure(figsize=(8, 5))
plt.plot(x_plot, y_plot, label="y(x)")
for b in region_boundaries[1:-1]:
    plt.axvline(b, color='gray', linestyle='--', alpha=0.6)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Multi-region BVP Solution")
plt.legend()
plt.grid(True)
plt.show()

How to Use This Template

Set region_boundaries — list of points where regions start/end.
Set region_k or your own parameters — one per region.
Modify fun(x, y) — define your ODE for each region.
Modify bc_multi — set your left/right boundary conditions
 and interface continuity rules.
Run — the solver will handle any number of regions.


Why This Works

Mesh Splitting: We explicitly split the mesh at each region boundary.
Continuity Enforcement: bc_multi enforces yyy and y′y'y′
continuity at each interface.
Flexible ODE: Inside fun, we detect which region each xxx
belongs to and apply the correct equation.
Scalable: Works for 2, 3, or more regions without changing the solver logic.


If you want, I can extend this so that each region can have
different ODE functions entirely (not just different parameters).
That would make it a fully modular multi-physics BVP solver.
Do you want me to prepare that version?
"""
# %% beam bending
def beam_bending_with_elastic_supports(bc,L=L, fun=None,
                               ec=None,gd=None,b=None,
                               point_loads=None,
                               q_func=None,
                               n_points=10):
    can_provide_total_load=False
    if fun is None:
        can_provide_total_load=True
        if q_func is None:
            q_func = lambda x: np.zeros_like(x)
        def q_total(x):
            q = q_func(x)
            if point_loads:
                for xp, P in point_loads:
                    q += P * (np.exp(-((x-xp)**2)/(2*1e-4))
                              / np.sqrt(2*np.pi*1e-4))
            return q
        if gd is None and b is None:
            def bending_fun(x, y):
                return np.vstack((y[1], y[2], y[3], q_total(x)/ec))
            fun=bending_fun
        else:
            def gbt_fun(x, y):
                y1, y2, y3, y4 = y
                dy1 = y2
                dy2 = y3
                dy3 = y4
                dy4 = (q_total(x) - gd * y3 - b * y1) / ec
                return np.vstack((dy1, dy2, dy3, dy4))
            fun=gbt_fun
    xs = np.linspace(0, L, n_points)
    ya = np.zeros((4,xs.size))
    sol=solve_bvp(fun, bc, xs, ya)
    if can_provide_total_load:
        total_load = np.trapezoid(q_total(xs), xs)
        if point_loads:
            total_load += sum(P for _, P in point_loads)
    else:
        total_load=None
        q_total=None
    return (sol,xs,q_total,total_load)

def run_example(bc,point_loads=None,
                gd=None,b=None,
                distributed_load=None,uc=None,L=L,fun=None,q=None,
                n_solve=10,n_plot=50):
    (res, xs, q_total, total_load) = (
        beam_bending_with_elastic_supports(
        bc, L=L, ec=E*Ixx,gd=gd,b=b,
        point_loads=point_loads,
        q_func=distributed_load,
        fun=fun,
        n_points=n_solve
    ))
    if q_total:
        q_plot=q_total
    else:
        q_plot=q
    if not res.success:
        raise Exception(f"solution failed, {res.message}")
    w_plot, theta_plot,M_plot, V_plot = [], [], [], []
    if n_solve == n_plot:
        px=xs
    else:
        px=np.linspace(0,L,n_plot)
    for x in px:
        y, dy, M, V = res.sol(x)
        w_plot.append(y)
        theta_plot.append(dy)
        M_plot.append(M)
        V_plot.append(V)

    # Deflection
    fig,axes=plt.subplot_mosaic(
    [
        ["q", None],
        ["w","theta"],
        ["M", "V"],
    ]
    ,num=f'solve bvp results for {uc}',clear=True)
    ax=axes["q"]
    ax.plot(xs, q_plot(xs), color="green", label="Load")
    if point_loads:
        for xp, P in point_loads:
            ax.axvline(x=xp, color="red", linestyle=":", alpha=0.7)
    ax.set_ylabel("q [N/m]")
    ax=axes["w"]
    ax.plot(xs, w_plot, label="Deflection")
    ax.set_ylabel("w [m]")
    if point_loads:
        for xp, P in point_loads:
            ax.axvline(x=xp, color="red", linestyle=":",
                        alpha=0.7, label="Point load"
                        if xp==point_loads[0][0] else "")
    ax=axes["theta"]
    ax.plot(xs, theta_plot, label="Rotation")
    ax.set_ylabel("θ [rad]")
    ax=axes["M"]
    ax.plot(xs, M_plot, label="Moment")
    if point_loads:
        for xp, P in point_loads:
            ax.axvline(x=xp, color="red", linestyle=":", alpha=0.7)
    ax.set_ylabel("M [Nm]")
    ax=axes["V"]
    ax.plot(xs, V_plot, label="Shear")
    if point_loads:
        for xp, P in point_loads:
            ax.axvline(x=xp, color="red", linestyle=":", alpha=0.7)
    ax.set_ylabel("V [N]")
    for ax in fig.axes:
        if ax.get_label()=='None':
            ax.remove()
        else:
            ax.legend()
            ax.set_xlabel("x [m]")
    fig.tight_layout()
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
uc='cf_dl'
res[uc]=run_example(bc,
                    distributed_load=distributed_load,uc=uc)
# %% Clamped beam with load at free end
def bc(ya,yb):
    return np.array([
        ya[0],# displacement(0)=0
        ya[1],# rotation'(0)=0
        yb[2],# moment''(L)=0
        yb[3] # shear'''(L)=0
    ])
point_loads = [(1*L, -L*A*rho*g)]
distributed_load = None
uc='cf_pl'
res[uc]=run_example(bc,point_loads=point_loads,uc=uc)
uc='2L_cf_pl'
point_loads = [(2*L, -2*L*A*rho*g)]
res[uc]=run_example(bc,point_loads=point_loads,uc=uc,L=2*L)
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
uc='ss_dl'
res[uc]=run_example(bc,point_loads,distributed_load,uc)
uc='2L_ss_dl'
res[uc]=run_example(bc,point_loads,distributed_load,uc=uc,L=2*L)
# %% Clamped beam with own weight and displacement support at x=L
def bc(ya,yb):
    return np.array([
        ya[0],# displacement(0)=0
        ya[1],# rotation'(0)=0
        yb[2],# moment''(L)=0
        yb[0] # displacement(L)=0
    ])
uc='cf_dl_s@L'
distributed_load = lambda x: -A*rho*g*np.ones_like(x)
res[uc]=run_example(bc,distributed_load=distributed_load,uc=uc)
# %% Clamped beam with own weight and stiff spring at x=0.8*L
# This is not working properly
def bc(ya,yb):
    return np.array([
        ya[0],# displacement(0)=0
        ya[1],# rotation'(0)=0
        yb[2],# moment''(L)=0
        yb[3] # shear'''(L)=0
    ])
ks=-A*rho*g*L/1e-6# Total load will cause at most 1e-6 displacement
def q(x):
    return -A*rho*g*np.ones_like(x)
def fun(x, y):
    mask = (x >= 0.799*L) & (x<=0.801*L)
    k = np.where(mask, ks, 0)
    y1, y2, y3, y4 = y
    dy1 = y2
    dy2 = y3
    dy3 = y4
    dy4 = (q(x) - k * y1) / (E*Ixx)
    return np.vstack((dy1, dy2, dy3, dy4))
uc='cf_dl_ss@0.8L'
res[uc]=run_example(bc,uc=uc,fun=fun,q=q)
# %% Clamped beam with own weight and displacement support at x=0.8*L
def bc(ya,yb):
    return np.array([
        ya[0],# displacement(0)=0
        ya[1],# rotation'(0)=0
        yb[2],# moment''(L)=0
        yb[0] # displacement(L)=0
    ])
uc='cf_dl_ds@L'
dl = lambda x: -A*rho*g*np.ones_like(x)
res[uc]=run_example(bc,distributed_load=dl,uc=uc)

