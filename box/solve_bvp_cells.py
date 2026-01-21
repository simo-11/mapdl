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
It = 10.513e-6
Iw = 1.64e-9
L = 0.15
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
def beam_bending_with_elastic_supports(L, EI, supports,
                               bc_left=None, bc_right=None,
                               point_loads=None, q_func=None,
                               n_points=20):
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
    
    def solve_segment(xa, xb, ya, yb, n_points=20):
        x = np.linspace(xa, xb, n_points)
        y_init = np.zeros((4, x.size))
    
        # Boundary condition function: ya_ and yb_ come from solve_bvp,
        # ya and yb are the target values you pass in
        def bc_segment(ya_, yb_, p=None):
            return np.array([
                ya_[0] - ya[0],   # deflection at left end
                ya_[1] - ya[1],   # slope at left end
                yb_[0] - yb[0],   # deflection at right end
                yb_[1] - yb[1]    # slope at right end
            ]) 
        sol = solve_bvp(bending_fun, bc_segment, x, y_init)
        return sol
    solutions = []
    ya = np.zeros(4)
    x_positions = [0] + [s[0] for s in supports] + [L]

    for i in range(len(x_positions)-1):
        xa, xb = x_positions[i], x_positions[i+1]
        yb_guess = np.zeros(4)
        sol = solve_segment(xa, xb, ya, yb_guess)
        solutions.append(sol)

        if i < len(supports):
            xm, k = supports[i]
            ym = sol.sol(xb)
            ym[2] = k/EI * ym[0]
            ya = ym

    sol_last = solutions[-1]
    ym = sol_last.sol(L)
    if bc_right:
        for cond, val in bc_right:
            if cond == "y": ym[0] = val
            elif cond == "dy": ym[1] = val
            elif cond == "d2y": ym[2] = val
            elif cond == "d3y": ym[3] = val

    def beam_response(x):
        for sol in solutions:
            if sol.x[0] <= x <= sol.x[-1]:
                vals = sol.sol(x)
                y, dy, d2y, d3y = vals
                M = EI * d2y
                V = EI * d3y
                return y, dy, M, V
        return None

    reactions = []
    for xm, k in supports:
        y, _, _, _ = beam_response(xm)
        R = k * y
        reactions.append((xm, R))

    xs = np.linspace(0, L, 500)
    total_load = np.trapezoid(q_total(xs), xs)
    if point_loads:
        total_load += sum(P for _, P in point_loads)
    total_reaction = sum(R for _, R in reactions)
    equilibrium_error = total_load + total_reaction

    return (beam_response, solutions, reactions, xs, q_total,
            total_load, equilibrium_error)
# %% beam ending example usage
L = 1.0
EI = 1.0
supports = [(0.3, 200.0), (0.6, 150.0)]
point_loads = [(0.5, -100.0)]
q_func = lambda x: np.ones_like(x)

(beam_response, solutions, reactions, xs, q_total, total_load,
 equilibrium_error) = (
    beam_bending_with_elastic_supports(
    L, EI, supports,
    bc_left=[("y",0),("dy",0)],
    bc_right=[("d2y",0),("d3y",0)],
    point_loads=point_loads,
    q_func=q_func
))

# Plot diagrams with markers and reaction annotations
x_plot = np.linspace(0, L, 200)
y_plot, M_plot, V_plot = [], [], []
for x in x_plot:
    y, dy, M, V = beam_response(x)
    y_plot.append(y)
    M_plot.append(M)
    V_plot.append(V)

plt.figure(figsize=(12,10))

# Deflection
plt.subplot(4,1,1)
plt.plot(x_plot, y_plot, label="Deflection")
for xm, _ in supports:
    plt.axvline(x=xm, color="blue", linestyle="--", 
                alpha=0.5, label="Support" if xm==supports[0][0] else "")
for xp, P in point_loads:
    plt.axvline(x=xp, color="red", linestyle=":", 
                alpha=0.7, label="Point load" if xp==point_loads[0][0] else "")
plt.ylabel("y(x)")
plt.title("Beam response diagrams")
plt.legend()

# Annotate reactions
for xm, R in reactions:
    plt.annotate(f"R={R:.1f}", xy=(xm,0), xytext=(xm,0.05),
                 arrowprops=dict(facecolor='black', shrink=0.05),
                 ha='center')

# Moment
plt.subplot(3,1,2)
plt.plot(x_plot, M_plot, label="Moment")
for xm, _ in supports:
    plt.axvline(x=xm, color="blue", linestyle="--", alpha=0.5)
for xp, P in point_loads:
    plt.axvline(x=xp, color="red", linestyle=":", alpha=0.7)
plt.ylabel("M(x)")

# Shear
plt.subplot(3,1,3)
plt.plot(x_plot, V_plot, label="Shear")
for xm, _ in supports:
    plt.axvline(x=xm, color="blue", linestyle="--", alpha=0.5)
for xp, P in point_loads:
    plt.axvline(x=xp, color="red", linestyle=":", alpha=0.7)
plt.ylabel("V(x)")
plt.xlabel("x")

# Continuous load visualization
plt.subplot(4,1,4)
plt.plot(xs, q_total(xs), color="green", label="Distributed load q(x)")
for xp, P in point_loads:
    plt.axvline(x=xp, color="red", linestyle=":", alpha=0.7)
plt.ylabel("q(x)")
plt.xlabel("x")
plt.legend()

plt.tight_layout()
plt.show()

print("Support reactions:")
for xm, R in reactions:
    print(f"Support at x={xm:.2f}, reaction={R:.2f}")
print(f"Equilibrium error (should be ~0): {equilibrium_error:.4e}")
