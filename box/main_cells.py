# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 10:48:02 2025

@author: simo nikula

Create, solve and report results using ansys models for cantilever BOX
"""
# %% commons
from concurrent.futures import ThreadPoolExecutor, TimeoutError
from ansys.mapdl import core as pymapdl
from ansys.dpf import core as dpf
from ansys.dpf import post
from ansys.dpf.core import operators#noqa
import matplotlib.pyplot as plt
import numpy as np
import types
import copy
import re
import enum
import psutil    
import time
class Model(enum.Enum):
    BEAM=1
    SOLID=2

def kill_mapdl_processes():
    for proc in psutil.process_iter(['pid', 'name']):
        if 'ansys' in proc.info['name'].lower():
            print(f"Killing PID {proc.info['pid']}, {proc.info['name']}")
            proc.kill()

def wait_for_shutdown(timeout=5):
    start = time.time()
    while time.time() - start < timeout:
        still_running = any('ansys' in p.info['name'].lower()
                            for p in psutil.process_iter(['name']))
        if not still_running:
            print("MAPDL shutdown confirmed.")
            return True
        time.sleep(0.5)
    print("Timeout reached — MAPDL may still be running.")
    return False


def stop_ansys():
    if 'mapdl' in globals():
        # Step 1: Try graceful shutdown
        try:
            mapdl.exit()
            del mapdl
        except Exception as e:
            print(f"Graceful exit failed: {e}")
            del mapdl
    # Step 2: Wait and verify
    if not wait_for_shutdown(timeout=5):
        # Step 3: Force kill if needed
        kill_mapdl_processes()  

"""
rotations are not supported in 0.14 (stable as of 2025-09)
so getting results using legacy methods
"""  
def pick_results(mapdl):
    mapdl.run("/solu")
    mapdl.antype("static")
    mapdl.nlgeom(key="on")
    mapdl.solve()
    mapdl.finish()
    mapdl.post1()
    mapdl.set(1)
    node_ids=mapdl.mesh.nnum
    node_coords=mapdl.mesh.nodes
    disp_data=mapdl.result.nodal_displacement(0)
    sns=types.SimpleNamespace(
        nnum=copy.deepcopy(node_ids),
        coords=copy.deepcopy(node_coords),
        displacement=copy.deepcopy(disp_data)
        )
    return sns

def get_sec_property(name):
    m=re.search(f'{name}\\s*=\\s*([-+\\d\\.E]+)',secdata)
    return float(m.group(1))


def check_mapdl(mapdl):
    return mapdl.is_alive  # Lightweight ping

if 'mapdl' in globals():
    try:
        with ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(check_mapdl, mapdl)#noqa
            try:
                is_alive = future.result(timeout=2)
                mapdl.clear()#noqa
                print(f"MAPDL is alive. Version: {mapdl.version}")#noqa
            except TimeoutError:
                print("MAPDL check timed out after 2 seconds.")
                del mapdl
            except Exception as e:
                print(f"MAPDL check failed: {e}")
                del mapdl
    except pymapdl.errors.MapdlExitedError:
        del mapdl
if not 'mapdl' in globals():
    mapdl=pymapdl.launch_mapdl()
    version=mapdl.version
    print(f"MAPDL lauched. Version: {version}")
# %% debug functions
# %% settings
models=(Model.BEAM,)
E=210E9
L=2
ndiv=20
master=0
rotation_in_degress=True
E=1E9
h=0.08
w=0.16
t=0.01
moment=1000
force=None
force_y=None
nu=0.3
G=E/(2*(1+nu))
sharp_corners=True
do_plots=False
if not do_plots:
    plt.close('all')
print(f"""BOX with h={h}, w={w}, t={t}, L={L} is active
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
mapdl.sectype(secid,"BEAM","HREC",'BOX',5)
secdata=mapdl.secdata(w,h,t,t,t,t)
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
mapdl.finish()
print("Beam model created")
# %% beam horizontal force
if force_y and Model.BEAM in models:
    mapdl.prep7()
    mapdl.fkdele('ALL','ALL')
    mapdl.fk(2,"FY",force_y)
    r_bhf=pick_results(mapdl)
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
    r_bt=pick_results(mapdl)
    print("Torsional load for beam processed")
# %% Solid model
if Model.SOLID in models:
    # https://ansyshelp.ansys.com/public/account/secured?returnurl=/////Views/Secured/corp/v242/en/ans_elem/Hlp_E_SOLID187.html
    mapdl.clear()
    mapdl.prep7()
    outer=mapdl.blc4(0,0,w,h,depth=L)
    inner=mapdl.blc4(t,t,w-2*t,h-2*t,depth=L)
    mapdl.vsbv(outer,inner,keep1='delete',keep2='delete')
    esize=np.pow(((2*w+2*h)*t*L)/1_000,1/3)
    if do_plots:
        mapdl.vplot(show_lines=True, 
                    line_width=5, 
                    show_bounds=True, 
                    cpos="iso")
    mapdl.et(1,"SOLID187")
    mapdl.mp("EX",1,E)
    mapdl.mp("PRXY",1,nu)
    mapdl.esize(esize)
    mapdl.vmesh('all')
    mapdl.nsel("S", "LOC", "Z", 0, 0)
    mapdl.d("ALL", "ALL", 0)
    mapdl.allsel()
    if do_plots:
        mapdl.eplot(plot_bc=True)
    print("Solid model created with "
          +f"{mapdl.mesh.n_elem} elements and {mapdl.mesh.n_node} nodes")
# %% torsion for solid using cerig
if Model.SOLID in models:
    mapdl.prep7()
    mapdl.fkdele('ALL','ALL')
    mapdl.cedele('all')
    mapdl.nsel("S", "LOC", "Z", L, L)       
    mapdl.nsel("R", "LOC", "X", w/2, w/2)       
    mapdl.nsel("R", "LOC", "Y", h/2, h/2) 
    mapdl.ndele('all')
    mapdl.allsel()
    mapdl.esel("S","ENAME","","MASS21")
    mapdl.edele('ALL')
    mapdl.allsel()
    master=np.max(mapdl.mesh.nnum)+1
    mapdl.n(master,w/2,h/2,L)
    mapdl.et(2,'MASS21')
    mapdl.type(2)
    mapdl.tshap('POINT')
    mapdl.r(1)
    mapdl.e(master)
    mapdl.nerr(nmerr=3,nmabt=1000_000)
    mapdl.nsel("S", "LOC", "Z", L, L)
    mapdl.cerig(master,'ALL',"UXYZ")
    mapdl.f(master,"MZ",moment)
    mapdl.allsel()
    r_st=pick_results()
    print("Torsional load for solid processed")
    if do_plots:
        model = dpf.Model(mapdl.result_file)
        disp = model.results.displacement()
        mesh=model.metadata.meshed_region
        fields_container = disp.outputs.fields_container()
        field = fields_container[0]
        mesh.plot(field,deform_by=disp, scale_factor=5.)
        sim=post.load_simulation(mapdl.result_file)
        displacement = sim.displacement()
        displacement.plot()
# %% plot results
def get_sorted_node_numbers(result):
    nnum=result.mesh.nnum
    nodes=result.mesh.nodes
    ss=nnum.tolist()
    ss.sort(key=lambda c:nodes[c-1][0])
    return ss

def plot_result(fig,ax,result,index,label):
    nnum =result.nnum
    disp = result.nodal_displacement  # disp shape: (n_nodes, 7)
    # Get coordinates
    coords = result.coords  # shape: (n_nodes, 3)
    # Extract X and ROTX (index 3), convert to degrees
    x_rotx = []
    for i, node in enumerate(nnum):
        x = coords[i][0]
        rotx_rad = disp[i][3]
        rotx_deg = rotx_rad * (180 / np.pi)
        x_rotx.append((x, rotx_deg))  
    # Sort by X
    x_rotx_sorted = sorted(x_rotx, key=lambda pair: pair[0])
    x_vals, rotx_vals = zip(*x_rotx_sorted)
    # Plot
    ax.plot(x_vals, rotx_vals, label=label)

def plot_solid_result(fig,ax,result,index,label):
    raise Exception("Not done")
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
    if rotation_in_degress:
        yv=180/np.pi*yv
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
    if rotation_in_degress:
        ax_t.set_ylabel(r'rotation [degrees]')
    else:
        ax_t.set_ylabel(r'rotation [radians]')
    if 'r_bt' in vars():
        plot_result(fig_t,ax_t,r_bt,3,'beam188')
    if 'r_st' in vars() and False:
        plot_solid_result(fig_t,ax_t,r_st,5,'solid187')
    add_analytical_torsion(ax_t,
                           get_sec_property('Torsion Constant'),
                           get_sec_property('Warping Constant'))
    ax_t.legend()
