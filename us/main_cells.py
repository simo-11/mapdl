# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 10:48:02 2025

@author: simo nikula

Create, solve and report results using ansys models for cantilever U-section
"""
# %% settings and commons initializations
from ansys.mapdl import core as pymapdl
import matplotlib.pyplot as plt
import numpy as np
h=0.1
w=0.05
t=0.004
E=210E9
nu=0.3
sharp_corners=True
do_plots=True
if not do_plots:
    plt.close('all')
L=2
ndiv=20
force=-1180
force_y=t/2
# moment for beam models must be set based on secdata
# and location of force
moment=None
no_moment_defined="moment is not defined for current geometry"
if h==0.1 and w == 0.05 and t == 0.004 and sharp_corners:
    moment=force*(0.015869+force_y)
else:
    print(no_moment_defined)
try:
    mapdl
except NameError:
    mapdl=pymapdl.Mapdl()
print(f"U with h={h}, w={w}, t={t} defined, do_plots={do_plots}")
# %% Beam model
mapdl.clear()
mapdl.prep7()
mapdl.et(1,"BEAM188")
mapdl.keyopt(1,1,1)#Warping included
#mapdl.keyopt(1,3,3)#Cubic shape functions
#mapdl.keyopt(1,4,2)#Combined shear stresses 
#mapdl.keyopt(1,7,2)#Output stresses at integration points
#mapdl.keyopt(1,9,3)#Output extrapolated stresses at nodes
mapdl.mp("EX",1,E)
mapdl.mp("PRXY",1,nu)
secid=1
mapdl.sectype(secid,"BEAM","CHAN","U",5)
mapdl.secoffset("CENT")
secdata=mapdl.secdata(w,w,h,t,t,t)
mapdl.k(1)
mapdl.k(2,L)
mapdl.lstr(1,2)
mapdl.lesize("ALL",ndiv=ndiv)
mapdl.lmesh("ALL")
mapdl.dk(kpoi=1,lab="ALL",value=0)
if do_plots:
    mapdl.nplot(vtk=True, nnum=True, cpos="xy",
            plot_bc=True,plot_bc_legend=True,
            bc_labels="mechanical",
            show_bounds=True, point_size=10)
print("Beam model created")
# %% beam horizontal force
mapdl.prep7()
mapdl.fkdele('ALL','ALL')
mapdl.fk(2,"FY",force)
mapdl.run("/solu")
mapdl.antype("static")
mapdl.solve()
r_bhf=mapdl.result
print("Horizontal load processed")
# %% beam vertical force
mapdl.prep7()
mapdl.fkdele('ALL','ALL')
mapdl.fk(2,"FZ",force)
mapdl.run("/solu")
mapdl.antype("static")
mapdl.solve()
r_bvf=mapdl.result
print("Vertical load processed")
# %% beam torsion
if moment is None:
    raise Exception(no_moment_defined)
mapdl.prep7()
mapdl.fkdele('ALL','ALL')
mapdl.fk(2,"MX",moment)
mapdl.run("/solu")
mapdl.antype("static")
mapdl.solve()
r_bt=mapdl.result
print("Torsional load processed")
# %% plot results
def get_sorted_node_numbers(result):
    nnum=result.mesh.nnum
    nodes=result.mesh.nodes
    ss=nnum.tolist()
    ss.sort(key=lambda c:nodes[c-1][0])
    return ss

def plot_result(fig,ax,result,index,label):
    sorted_node_numbers=get_sorted_node_numbers(result)
    size=result.mesh.nnum.size
    xvs=result.mesh.nodes[:,0]
    nd=result.nodal_displacement(0)[1][:,index]
    xv=np.zeros(size).tolist()
    yv=np.zeros(size).tolist()
    i=0
    for nn in sorted_node_numbers:
        ni=nn-1
        xv[i]=xvs[ni]
        yv[i]=abs(nd[ni])
        i=i+1
    ax.plot(xv,yv,label=label)
    ax.legend()
fig_hf, ax_hf = plt.subplots(num='horizontal force',clear=True)
ax_hf.set_xlabel(r'x-coordinate [m]')
ax_hf.set_ylabel(r'horizontal displacement [m]')
plot_result(fig_hf,ax_hf,r_bhf,1,'beam')
fig_vf, ax_vf = plt.subplots(num='vertical force',clear=True)
ax_vf.set_xlabel(r'x-coordinate [m]')
ax_vf.set_ylabel(r'vertical displacement [m]')
plot_result(fig_vf,ax_vf,r_bvf,2,'beam')
fig_t, ax_t = plt.subplots(num='torsion',clear=True)
ax_t.set_ylabel(r'x-coordinate [m]')
ax_t.set_ylabel(r'rotation [radians]')
plot_result(fig_t,ax_t,r_bt,3,'beam')
