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
import types
import copy
import re
def pick_results():
    mapdl.run("/solu")
    mapdl.antype("static")
    mapdl.solve()
    r=mapdl.result
    return types.SimpleNamespace(
        mesh=copy.deepcopy(r.mesh), 
        nodal_displacement=copy.deepcopy(r.nodal_displacement(0)))

def get_sec_property(name):
    m=re.search(f'{name}\\s*=\\s*([-+\\d\\.E]+)',secdata)
    return float(m.group(1))

h=0.1
w=0.05
t=0.004
E=210E9
nu=0.3
G=E/(2*(1+nu))
sharp_corners=True
do_plots=False
if not do_plots:
    plt.close('all')
L=2
ndiv=20
force=-1180
force_y=t/2
try:
    mapdl
    try:
        mapdl.clear()
    except pymapdl.errors.MapdlExitedError:
        del mapdl        
except NameError:
    try:
        mapdl=pymapdl.Mapdl(timeout=5)
    except pymapdl.errors.MapdlConnectionError:
        mapdl=pymapdl.launch_mapdl();
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
#mapdl.secoffset("SHRC")#CENT,SHRC,ORIGIN,USER
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
r_bhf=pick_results()
print("Horizontal load processed")
# %% beam vertical force
# add moment to cancel effect of torsion
moment=-force*(
    get_sec_property('Centroid Y')-get_sec_property('Shear Center Y'))
mapdl.prep7()
mapdl.fkdele('ALL','ALL')
mapdl.fk(2,"FZ",force)
mapdl.fk(2,"MX",moment)
r_bvf=pick_results()
print("Vertical load processed")
# %% beam torsion if vertical load applied at force_y
moment=force*(force_y-get_sec_property('Shear Center Y'))
mapdl.prep7()
mapdl.fkdele('ALL','ALL')
mapdl.fk(2,"MX",moment)
r_bt=pick_results()
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
    nd=result.nodal_displacement[1][:,index]
    xv=np.zeros(size).tolist()
    yv=np.zeros(size).tolist()
    i=0
    for nn in sorted_node_numbers:
        ni=nn-1
        xv[i]=xvs[ni]
        yv[i]=nd[ni]
        i=i+1
    ax.plot(xv,yv,label=label)

def add_analytical_bending(ax,I):
    xv=np.linspace(0,L)
    yv=force*np.pow(L,3)/(6*E*I)*(2-3*(L-xv)/L+np.pow((L-xv)/L,3))
    ax.plot(xv,yv,label='analytical')

def kc(It,Iw):
    return np.sqrt((G*It)/(E*Iw))

def theta(T,It,Iw,L,x):
    k=kc(It,Iw)
    c0=T/(k*G*It)
    if np.tanh(k*L)==1:
        y=c0*(np.exp(-k*x)-1+k*x)
    else:
        y=c0*((np.tanh(k*L)*(np.cosh(k*x)-1))-np.sinh(k*x)+k*x)
    return y

def add_analytical_torsion(ax,It,Iw):
    xv=np.linspace(0,L)
    yv=theta(moment,It,Iw,L,xv)
    ax.plot(xv,yv,label='analytical')
    
fig_hf, ax_hf = plt.subplots(num='horizontal force',clear=True)
ax_hf.set_xlabel(r'x-coordinate [m]')
ax_hf.set_ylabel(r'horizontal displacement [m]')
plot_result(fig_hf,ax_hf,r_bhf,1,'beam188')
add_analytical_bending(ax_hf,get_sec_property('Izz'))
ax_hf.legend()
fig_vf, ax_vf = plt.subplots(num='vertical force',clear=True)
ax_vf.set_xlabel(r'x-coordinate [m]')
ax_vf.set_ylabel(r'vertical displacement [m]')
plot_result(fig_vf,ax_vf,r_bvf,2,'beam188')
add_analytical_bending(ax_vf,get_sec_property('Iyy'))
ax_vf.legend()
fig_t, ax_t = plt.subplots(num='torsion',clear=True)
ax_t.set_ylabel(r'x-coordinate [m]')
ax_t.set_ylabel(r'rotation [radians]')
plot_result(fig_t,ax_t,r_bt,3,'beam188')
add_analytical_torsion(ax_t,
                       get_sec_property('Torsion Constant'),
                       get_sec_property('Warping Constant'))
ax_t.legend()
