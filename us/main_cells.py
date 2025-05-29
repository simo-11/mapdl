# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 10:48:02 2025

@author: simo nikula

Create, solve and report results using ansys models for cantilever 
 * U-section
 * BOX
"""
# %% commons
from ansys.mapdl import core as pymapdl
import matplotlib.pyplot as plt
import numpy as np
import types
import copy
import re
import enum

class Section(enum.Enum):
    U=1
    BOX=2

class Model(enum.Enum):
    BEAM=1
    SOLID=2    
    
def pick_results():
    mapdl.run("/solu")
    mapdl.antype("static")
    mapdl.nlgeom(key="on")
    mapdl.solve()
    r=mapdl.result
    return types.SimpleNamespace(
        mesh=copy.deepcopy(r.mesh), 
        nodal_displacement=copy.deepcopy(r.nodal_displacement(0)))

def get_sec_property(name):
    m=re.search(f'{name}\\s*=\\s*([-+\\d\\.E]+)',secdata)
    return float(m.group(1))

if 'mapdl' in vars():
    try:
        mapdl.clear()#noqa
    except pymapdl.errors.MapdlExitedError:
        del mapdl        
if not 'mapdl' in vars():
    try:
        mapdl=pymapdl.Mapdl(timeout=5)
    except pymapdl.errors.MapdlConnectionError:
        mapdl=pymapdl.launch_mapdl()

# %% settings
section=Section.BOX
models=(Model.SOLID,)
E=210E9
L=2
ndiv=20
match section.name:
    case "U":
        h=0.1
        w=0.05
        t=0.004
        force=-1180
        force_y=t/2
        moment=None
    case "BOX":
        E=1E9
        h=0.08
        w=0.16
        t=0.01
        moment=1000
        force=None
    case _:
        raise Exception(f'parameters for {section.name} are not defined')    
nu=0.3
G=E/(2*(1+nu))
sharp_corners=True
do_plots=True
if not do_plots:
    plt.close('all')
print(f"""{section.name} with h={h}, w={w}, t={t} is active
E={E:.3G}, nu={nu:.3G}
models={models}
do_plots={do_plots}
pymapdl_version={mapdl.info._get_pymapdl_version()}""")
# %% beam model
# secdata is needed for analytical solution
mapdl.clear()
mapdl.prep7()
mapdl.et(1,"BEAM188")
mapdl.keyopt(1,1,1)
# 1=Warping included
#mapdl.keyopt(1,2,0)
# Cross-section scaling for nlgeom
#mapdl.keyopt(1,3,3)
# Cubic shape functions, 0=Linear, 1=Quadratic
#mapdl.keyopt(1,4,0)
# 0=torsion related, 1=transverese, 2=Combined shear stresses 
#mapdl.keyopt(1,5,0)
# 0=3D, 1=XY
#mapdl.keyopt(1,6,0)
# 0=Output section forces/moments and strains/curvatures at 
#   integration points along the length (default)
#mapdl.keyopt(1,7,0)
# 0=None, 1=Max and min, 2=Output stresses at each section point
#mapdl.keyopt(1,9,0)
# 0=None
# 1=Maximum and minimum stresses/strains
# 2=1 plus stresses and strains along the exterior 
#   boundary of the cross-section
# 3=1 plus stresses and strains at all section nodes
#mapdl.keyopt(1,11,0)
# 0=Automatically determine if preintegrated section properties can be used
# 1=Use numerical integration of section
#mapdl.keyopt(1,12,0)
# 0=Linear tapered section analysis; cross-section properties 
#   are evaluated at each Gauss point
# 1=Average cross-section analysis
#mapdl.keyopt(1,15,0)
# 0=Store averaged results at each section corner node
# 1=Store non-averaged results at each section integration point.
mapdl.mp("EX",1,E)
mapdl.mp("PRXY",1,nu)
secid=1
match section.name:
    case "U":
        mapdl.sectype(secid,"BEAM","CHAN",section.name,5)
        #mapdl.secoffset("SHRC")#CENT,SHRC,ORIGIN,USER
        secdata=mapdl.secdata(w,w,h,t,t,t)
    case "BOX":
        mapdl.sectype(secid,"BEAM","HREC",section.name,5)
        secdata=mapdl.secdata(w,h,t,t,t,t)
    case _:
        raise Exception(f'parameters for {section.name} are not defined')    
mapdl.k(1)
mapdl.k(2,L)
mapdl.lstr(1,2)
mapdl.lesize("ALL",ndiv=ndiv)
mapdl.lmesh("ALL")
mapdl.dk(kpoi=1,lab="ALL",value=0)
if do_plots:
    mapdl.nplot(nnum=True, cpos="xy",
            plot_bc=True,plot_bc_legend=True,
            bc_labels="mechanical",
            show_bounds=True, point_size=10)
print("Beam model created")
# %% beam horizontal force
if force and Model.BEAM in models:
    mapdl.prep7()
    mapdl.fkdele('ALL','ALL')
    mapdl.fk(2,"FY",force)
    r_bhf=pick_results()
    print("Horizontal load for beam processed")
# %% beam vertical force
# add moment to cancel effect of torsion
if force and Model.BEAM in models:
    moment=-force*(
        get_sec_property('Centroid Y')-get_sec_property('Shear Center Y'))
    mapdl.prep7()
    mapdl.fkdele('ALL','ALL')
    mapdl.fk(2,"FZ",force)
    mapdl.fk(2,"MX",moment)
    r_bvf=pick_results()
    print("Vertical load for beam processed")
# %% beam torsion 
# if moment is not given, uses vertical load applied at force_y
if Model.BEAM in models:
    if not moment:
        moment=force*(force_y-get_sec_property('Shear Center Y'))
    mapdl.prep7()
    mapdl.fkdele('ALL','ALL')
    mapdl.fk(2,"MX",moment)
    r_bt=pick_results()
    print("Torsional load for beam processed")
# %% Solid model
if Model.SOLID in models:
    # https://ansyshelp.ansys.com/public/account/secured?returnurl=/////Views/Secured/corp/v242/en/ans_elem/Hlp_E_SOLID187.html
    mapdl.clear()
    mapdl.prep7()
    outer=mapdl.blc4(0,0,w,h,depth=L)
    inner=mapdl.blc4(t,t,w-2*t,h-2*t,depth=L)
    mapdl.vsbv(outer,inner,keep1='delete',keep2='delete')
    if do_plots:
        mapdl.vplot(show_lines=True, 
                    line_width=5, 
                    show_bounds=True, 
                    cpos="iso")
    mapdl.et(1,"SOLID187")
    mapdl.mp("EX",1,E)
    mapdl.mp("PRXY",1,nu)
    mapdl.esize(5*t)
    mapdl.vmesh('all')
    mapdl.nsel("S", "LOC", "Z", 0, 0)
    mapdl.d("ALL", "ALL", 0)
    mapdl.allsel()
    if do_plots:
        mapdl.eplot()
        mapdl.nplot(plot_bc=True)
    print("Solid model created with "
          +f"{mapdl.mesh.n_elem} elements and {mapdl.mesh.n_node} nodes")
# %% torsion for solid 
if Model.SOLID in models:
    print("Torsional load for solid processed")
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

if force:    
    fig_hf, ax_hf = plt.subplots(num='horizontal force',clear=True)
    ax_hf.set_xlabel(r'x-coordinate [m]')
    ax_hf.set_ylabel(r'horizontal displacement [m]')
    if 'r_bhf' in vars():
        plot_result(fig_hf,ax_hf,r_bhf,1,'beam188')
    add_analytical_bending(ax_hf,get_sec_property('Izz'))
    ax_hf.legend()
if force:    
    fig_vf, ax_vf = plt.subplots(num='vertical force',clear=True)
    ax_vf.set_xlabel(r'x-coordinate [m]')
    ax_vf.set_ylabel(r'vertical displacement [m]')
    if 'r_bvf' in vars():
        plot_result(fig_vf,ax_vf,r_bvf,2,'beam188')
    add_analytical_bending(ax_vf,get_sec_property('Iyy'))
    ax_vf.legend()
if moment:
    fig_t, ax_t = plt.subplots(num='torsion',clear=True)
    ax_t.set_ylabel(r'x-coordinate [m]')
    ax_t.set_ylabel(r'rotation [radians]')
    if 'r_bt' in vars():
        plot_result(fig_t,ax_t,r_bt,3,'beam188')
    add_analytical_torsion(ax_t,
                           get_sec_property('Torsion Constant'),
                           get_sec_property('Warping Constant'))
    ax_t.legend()
