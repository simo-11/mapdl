# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 10:48:02 2025

@author: simo nikula

Create, solve and report results using ansys models for cantilever BOX

To reuse same ansys mapdl instance it is used as global variable
and Spyder should be configured to Run in console's namespace
"""
#%% commons
from concurrent.futures import ThreadPoolExecutor, TimeoutError
from ansys.mapdl import core as pymapdl
from ansys.dpf import core as dpf
from ansys.dpf import post
from ansys.dpf.core import operators#noqa
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
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
    
def bm(keyopt1=0,keyopt3=0):
    mapdl.clear()
    mapdl.prep7()
    mapdl.et(1,"BEAM188")
    # https://ansyshelp.ansys.com/public/account/secured?returnurl=/Views/Secured/corp/v252/en/ans_elem/Hlp_E_BEAM188.html?q=beam188
    # Keyopt 1:
    # 0=No warping (default)
    # 1=Warping included
    #   Causes numerical issues if torsion exceeds about 2π 
    # 2=Remove warping for closed sections
    mapdl.keyopt(1,1,keyopt1)
    # Keyopt 3: 
    # shape functions along then length
    # 0=Linear (default), migitates numerical issues partly
    # 2=Quadratic 
    # 3=Qubic
    #
    mapdl.keyopt(1,3,keyopt3)
    mapdl.mp("EX",1,E)
    mapdl.mp("PRXY",1,nu)
    secid=1
    mapdl.sectype(secid,"BEAM","HREC",'BOX',5)
    secdata=mapdl.secdata(w,h,t,t,t,t)
    mapdl.k(1)
    mapdl.k(2,L)
    mapdl.lstr(1,2)
    mapdl.lesize("ALL",ndiv=ndiv,space=-40)
    mapdl.lmesh("ALL")
    mapdl.dk(kpoi=1,lab="ALL",value=0)
    if do_plots:
        mapdl.nplot(nnum=True, cpos="xy",
                plot_bc=True,plot_bc_legend=True,
                bc_labels="mechanical",
                show_bounds=True, point_size=10)
    return secdata

"""
rotations are not supported in 0.14 (stable as of 2025-09)
so getting results using legacy methods from 
uncompressed results file requested using /fcomp,rst,0
"""  
def pick_results(mapdl,nlgeom=False,file=None):
    if file!=None:
        mapdl.finish()
        mapdl.title(file)
        mapdl.filname(fname=file)
    mapdl.run("/solu")
    mapdl.run("outres,all,all")
    mapdl.antype("static")
    if nlgeom:
        mapdl.nlgeom(key="on")
    else:
        mapdl.nlgeom(key="off")
    solve_txt=mapdl.solve()
    if nlgeom:
        sol=mapdl.solution
        if not sol.converged:
            raise Exception("Solution failed, "
                    f"converged={sol.converged},"
                   f" moment_cnv={sol.moment_cnv:.5g},"
                   f" n_eqit={sol.n_eqit:.0f}")
    else:
        sol=None
    mapdl.finish()
    mapdl.post1()
    mapdl.set(1,1)
    node_coords=mapdl.mesh.nodes
    node_ids,disp_data=mapdl.result.nodal_displacement(0)
    sns=types.SimpleNamespace(
        file=file,
        nnum=copy.deepcopy(node_ids),
        coords=copy.deepcopy(node_coords),
        displacement=copy.deepcopy(disp_data),
        result_file=mapdl.result_file,
        solve_txt=solve_txt,
        sol=sol
        )
    return sns

def dpf1(sns, scale_factor=None):
    if scale_factor==None:
        scale_factor=1.
    model = dpf.Model(sns.result_file)
    disp = model.results.displacement()
    mesh=model.metadata.meshed_region
    fields_container = disp.outputs.fields_container()
    field = fields_container[0]
    mesh.plot(field,deform_by=disp, scale_factor=scale_factor)
    
def dpf2(sns):
    sim=post.load_simulation(sns.result_file)
    displacement = sim.displacement()
    displacement.plot()

def get_sec_property(secdata,name):
    m=re.search(f'{name}\\s*=\\s*([-+\\d\\.E]+)',secdata)
    return float(m.group(1))

def check_mapdl(mapdl):
    return mapdl.is_alive

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
        global mapdl
        # Step 1: Try graceful shutdown
        try:
            mapdl.exit()
        except Exception as e:
            print(f"Graceful exit failed: {e}")
    # Step 2: Wait and verify
    if not wait_for_shutdown(timeout=5):
        # Step 3: Force kill if needed
        kill_mapdl_processes()
    if 'mapdl' in globals():
        del globals()['mapdl']

def check_global_mapdl():
    if not 'mapdl' in globals():
        print("global mapdl is not defined, ansys process will be started")
        return
    global mapdl
    timeout=2
    with ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(check_mapdl, mapdl)#noqa
        try:
            is_alive = future.result(timeout=timeout)
            if is_alive:
                mapdl.clear()#noqa
                print(f"MAPDL is alive. Version: {mapdl.version}")#noqa
                return
        except TimeoutError:
            print(f"MAPDL check timed out after {timeout} seconds.")
        except Exception as e:
            print(f"MAPDL check failed: {e}")
        except pymapdl.errors.MapdlExitedError as e:
                print(f"MAPDL check exited: {e}")
    stop_ansys()

check_global_mapdl()

if not 'mapdl' in globals():
    try:
        mapdl=pymapdl.launch_mapdl(run_location="local",
                                   override=True,
                                   cleanup_on_exit=False)
    except Exception as e:
        print(f"Launch failed: {e}")
        stop_ansys()
        mapdl=pymapdl.launch_mapdl(run_location="local",
                                   override=True,
                                   cleanup_on_exit=False)
    mapdl.units("mks")
    mapdl.run("/FCOMP,RST,0")
    version=mapdl.version
    print(f"MAPDL lauched. Version: {version}")
models=(Model.BEAM,)
L=2
ndiv=10
master=0
rotation_in_degress=True
E=1E9
h=0.08
w=0.16
t=0.01
moment=1000
force=0#1000
force_y=0#1000
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
pymapdl_version={mapdl.info._get_pymapdl_version()}
moment={moment}, force={force}, force_y={force_y}
""")
#%% debug functions
#%% secdata
mapdl.clear()
mapdl.prep7()
for refinekey in range(0,6):
    secid=refinekey
    start=time.time()
    mapdl.sectype(secid,"BEAM","HREC",f'BOX{refinekey}',refinekey)
    sd=mapdl.secdata(w,h,t,t,t,t)
    end=time.time()
    It=get_sec_property(sd,'Torsion Constant')
    Iw=get_sec_property(sd,'Warping Constant')
    print("{0} & {1:.4g} & {2:.4g} & {3:.2g} \\\\".
          format(refinekey,
      1e6*It,
      1e9*Iw,
      (end-start)
      ))
et=mapdl.et(ename='plane183')
for ec_in_h in range(10,40,5):
    mapdl.asel('all')
    mapdl.aclear('all')
    mapdl.adele('all')
    anum0 = mapdl.blc4(0, 0, w, h)
    anum1 = mapdl.blc4(t, t, w-2*t, h-2*t)
    aout = mapdl.asba(anum0, anum1)
    name=f'BOX{ec_in_h}'
    fname=f'mesh_{name}'
    secid=secid+1
    start=time.time()
    lesize=h/ec_in_h
    mapdl.lesize('all',lesize)
    mapdl.amesh('all')
    nnodes = int(mapdl.get_value("NODES",0, "COUNT"))
    mapdl.secwrite(fname=fname,elem_type=et)
    mapdl.sectype(secid,"BEAM","MESH",name)
    mapdl.secread(fname=fname,option='mesh')
    end=time.time()
    It=mapdl.get_value('secp',secid,'prop','tors')
    Iw=mapdl.get_value('secp',secid,'prop','warp')
    print("{0} & {1:.4g} & {2:.4g} & {3:.2g} \\\\".
          format(nnodes,
      1e6*It,
      1e9*Iw,
      (end-start)
      ))    
#%% beam horizontal force
if force_y and Model.BEAM in models:
    bm()
    mapdl.fk(2,"FY",force_y)
    r_bhf=pick_results(mapdl,file='bhf')
    r_bhf_nl=pick_results(mapdl,True,file='bhf_nl')
    print("Horizontal load for beam processed")
#%% beam vertical force
# add moment to cancel effect of torsion
if force and Model.BEAM in models:
    secdata=bm()
    _moment=-force*(
        get_sec_property(secdata,'Centroid Y')
        -get_sec_property(secdata,'Shear Center Y'))
    mapdl.fk(2,"FZ",force)
    mapdl.fk(2,"MX",_moment)
    r_bvf=pick_results(mapdl,file='bvf')
    r_bvf_nl=pick_results(mapdl,True,file='bvf_nl')
    print("Vertical load for beam processed")
#%% beam torsion bt
if Model.BEAM in models:    
    _moment=moment
    bm(keyopt1=1)
    mapdl.dk(kpoi=2,lab="UY",lab2="UZ",value=0)    
    mapdl.fk(2,"MX",_moment)
    print(f"Torsional load of {_moment} for beam (bt) is processed")
    r_bt=pick_results(mapdl,file='bt')
    r_bt.rfe=r_bt.displacement[1][3]*180/np.pi
    print(f"Rotation for bt with nlgeom=off (bt) {r_bt.rfe:.4g}°")
    r_bt_nl=pick_results(mapdl,True,file='bt_nl')
    r_bt_nl.rfe=r_bt_nl.displacement[1][3]*180/np.pi
    print(f"Rotation for bt with nlgeom=on (bt_nl) {r_bt_nl.rfe:.4g}°")
#%% bt1 
# warping constrained at loaded end
if Model.BEAM in models:
    _moment=moment
    bm(keyopt1=1)
    mapdl.dk(kpoi=2,lab="UY",lab2="UZ",lab3="WARP",value=0)
    mapdl.fk(2,"MX",_moment)
    print(f"Torsional load of {_moment} for beam (bt1) is processed")
    r_bt1=pick_results(mapdl,file='bt1')
    r_bt1.rfe=r_bt1.displacement[1][3]*180/np.pi
    print(f"Rotation for bt1 with nlgeom=off (bt1) {r_bt1.rfe:.4g}°")
    r_bt1_nl=pick_results(mapdl,True,file='bt1_nl')
    r_bt1_nl.rfe=r_bt1_nl.displacement[1][3]*180/np.pi
    print(("Rotation for bt1 with nlgeom=on (bt1_nl)"
           f" {r_bt1_nl.rfe:.4g}°"))
#%% bt2 
# warping and axial displacement constrained at loaded end
uc='bt2'
if Model.BEAM in models:
    _moment=moment
    bm(keyopt1=1)
    mapdl.dk(kpoi=2,lab="UY",lab2="UZ",lab3="WARP",lab4="UX",value=0)  
    mapdl.fk(2,"MX",_moment)
    print(f"Torsional moment of {_moment} for beam ({uc}) is processed")
    r_bt2=pick_results(mapdl,file=f'{uc}')
    r_bt2.rfe=r_bt2.displacement[1][3]*180/np.pi
    print(f"Rotation for {uc} with nlgeom=off ({uc}) {r_bt2.rfe:.4g}°")
    r_bt2_nl=pick_results(mapdl,True,file=f'{uc}_nl')
    r_bt2_nl.rfe=r_bt2_nl.displacement[1][3]*180/np.pi
    print((f"Rotation for {uc} with nlgeom=on ({uc}_nl)"
           f" {r_bt2_nl.rfe:.4g}°"))
#%% bt3 
# warping excluded
uc='bt3'
if Model.BEAM in models:
    _moment=moment
    bm(keyopt1=0)
    mapdl.dk(kpoi=2,lab="UY",lab2="UZ",lab3="WARP",lab4="UX",value=0)  
    mapdl.fk(2,"MX",_moment)
    print(f"Torsional moment of {_moment} for beam ({uc}) is processed")
    r_bt2=pick_results(mapdl,file=f'{uc}')
    r_bt2.rfe=r_bt2.displacement[1][3]*180/np.pi
    print(f"Rotation for {uc} with nlgeom=off ({uc}) {r_bt2.rfe:.4g}°")
    r_bt2_nl=pick_results(mapdl,True,file=f'{uc}_nl')
    r_bt2_nl.rfe=r_bt2_nl.displacement[1][3]*180/np.pi
    print((f"Rotation for {uc} with nlgeom=on ({uc}_nl)"
           f" {r_bt2_nl.rfe:.4g}°"))
#%% btol
# displacement constraint at loaded end and keyopts varied
# if warping is included solution fails for rot_vals n*2*pi
# 
rot_vals = np.concatenate((np.linspace(0., 1.5 * np.pi, 5),
           np.linspace(1.65*np.pi, 2 * np.pi, 5)))
max_rot=max(rot_vals)
xmax=180/np.pi*max_rot
if not 'secdata' in globals():
    secdata=bm()
ymax=1.2*get_sec_property(secdata,'Torsion Constant')*G*max_rot/L
fig=plt.figure(num='btol')
fig.clear()
ax=fig.gca()
ax.clear()
ax.set_xlim(0, xmax)
ax.set_ylim(0, ymax)
ax2 = ax.twinx()
ax2.clear()
ax2.set_ylim(-ymax/2, 6*ymax)
plt.xlabel('Rotation [degrees]')
ax.set_ylabel('Reaction Moment [Nm]')
ax2.set_ylabel('Reaction Force [N]')
ax2.yaxis.set_label_position("right")
plt.title('Reaction vs. Rotation using Ansys Beams with NLGEOM)')
plt.grid(True)
plt.tight_layout()
for uc in range(1,3): # use upper limit of 4 or 5 to see failure
    match uc:
        case 1|2:
            keyopt1=0
            keyopt3=0
            marker=MarkerStyle((3,2,0))
            color='blue'
            color2='cyan'
            linestyle=(0,(2*uc,2*uc))
        case _ if uc<5:
            keyopt1=1
            keyopt3=uc-2
            marker=MarkerStyle((3,1,0))
            color='red'
            color2='orange'
            linestyle=(0,(2*uc,2*uc))
        case _:raise Exception(f'Settings for uc {uc} are missing')
    label=f"uc{uc}"        
    line,=ax.plot([],[],label=label,
         marker=marker, linestyle=linestyle, color=color)
    line2,=ax2.plot([],[],label=label,
         marker=marker, linestyle=linestyle, color=color2)
    ax.legend(loc='upper left')
    ax2.legend(loc='upper center')
    bm(keyopt1=keyopt1,keyopt3=keyopt3)
    if uc % 2 == 0:
        mapdl.dk(kpoi=2,lab="UY",lab2="UZ",lab3="WARP",lab4="UX",value=0) 
    else:
        mapdl.dk(kpoi=2,lab="UY",lab2="UZ",value=0) 
    print(f"btol {label} is processed")
    mapdl.finish()
    mapdl.title(label)
    mapdl.filname(fname=label)
    mapdl.run("/SOLU")
    mapdl.antype("STATIC")
    mapdl.nlgeom("ON")
    mapdl.nropt("FULL")
    mapdl.kbc(0) 
    mapdl.autots("OFF")
    mapdl.outres("ALL", "ALL")
    mapdl.autots()
    rotation_deg = [0]
    reaction_rx = [0]
    reaction_fx = [0]
    try:
        for i, rot in enumerate(rot_vals, start=1):
            mapdl.time(i)
            mapdl.d(2, "ROTX", rot)
            rv=f"{rot*180/np.pi:4.3g}°"
            print(f"Target rotation for step {i:2d} is {rot:4.3g} i.e. {rv}",
                  end=", ")      
            solve_text=mapdl.solve()
            sol=mapdl.solution
            if not sol.converged:
                print(f"converged={sol.converged},"
                       f" moment_cnv={sol.moment_cnv:.5g},"
                       f" n_eqit={sol.n_eqit:.0f}")
                break
            mx = np.abs(mapdl.get_value("NODE", 1, "RF", "mx"))
            fx = np.abs(mapdl.get_value("NODE", 1, "RF", "fx"))
            arot = mapdl.post_processing.nodal_rotation('x')[1]            
            ra = arot*180/np.pi
            print(f"rf/mx={mx:.5g}, rf/fx={fx:.5g}"
                   f" moment_cnv={sol.moment_cnv:.5g},"
                   f" n_eqit={sol.n_eqit:.0f}")
            reaction_rx.append(mx)
            reaction_fx.append(fx)
            rotation_deg.append(ra)
            line.set_data(rotation_deg, reaction_rx)
            line2.set_data(rotation_deg, reaction_fx)
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.pause(0.05)    
    except Exception as e:
        print(f"Solve failed for {rv} due to {e}")        
#%% qtplot
from pyvistaqt import BackgroundPlotter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as plt_colors
import pathlib

def qtplot(src,scale=None,scalar_component='UX',node_labels=False):
    # Load result file
    if isinstance(src,str):
        fn=f"local/{src}.rst"
        file=src
    elif (isinstance(src,types.SimpleNamespace) and 
        isinstance(src.result_file,pathlib.PurePath)):
        fn=src.result_file
        file=src.file
    else:
        raise TypeError(
            ("Expected src to be str "
             "or namespace having result_file "
             f"but got {type(src).__name__}"))
    model = dpf.Model(fn)
    # Get mesh and displacement field
    mesh = model.metadata.meshed_region
    disp_fc = model.results.displacement().eval()
    disp = disp_fc[0]  # Field
    # Convert mesh to PyVista format
    grid = mesh.grid
    # Build node ID → point index mapping
    node_ids = mesh.nodes.scoping.ids  # DPF node IDs
    point_id_map = {node_id: i for i, node_id in enumerate(node_ids)}
    # Create displacement array aligned with PyVista point order
    n_points = grid.n_points
    disp_array = np.zeros((n_points, 3))
    if node_labels:
        labels = np.empty(n_points,dtype='<U30')
    # Use disp.scoping.ids to get node IDs for each displacement vector
    for i, node_id in enumerate(disp.scoping.ids):
        if node_id in point_id_map:
            idx = point_id_map[node_id]
            disp_array[idx] = disp.data[i]
            if node_labels:
                labels[idx]=(f"{node_id} ("
                         f"{grid.points[idx][0]+disp_array[idx][0]:.4f}, "
                         f"{grid.points[idx][1]+disp_array[idx][1]:.4f}, "
                         f"{grid.points[idx][2]+disp_array[idx][2]:.4f})"
                         )
    # Apply scaled displacement to mesh
    if scale==None:
        # Compute bounding box dimensions
        bounds = grid.bounds
        x_range = bounds[1] - bounds[0]
        y_range = bounds[3] - bounds[2]
        z_range = bounds[5] - bounds[4]
        max_dim = max(x_range, y_range, z_range)        
        # Compute displacement magnitude
        disp_magnitude = np.linalg.norm(disp_array, axis=1)
        max_disp = disp_magnitude.max()
        # Autoscale: max displacement = 10% of largest model dimension
        target_disp = 0.1 * max_dim
        scale = target_disp / max_disp if max_disp > 0 else 1.0
    deformed_grid = grid.copy()
    deformed_grid.points = grid.points + disp_array * scale
    # Choose scalar component (e.g. UX)
    scalars = disp_array[:, 0]  # UX
    # Create a diverging colormap centered at zero
    cmap = plt.get_cmap("coolwarm")  # or "seismic", "RdBu", "PiYG", etc.
    vmin=scalars.min()
    vmax=scalars.max()
    if vmin<0 and vmax>0:
        # Normalize so that zero is white
        norm = plt_colors.TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    else:
        norm = plt_colors.Normalize(vmin=vmin, vmax=vmax)
    colors = cmap(norm(scalars))[:, :3]  # Drop alpha channel
    # Launch interactive non-blocking window
    plotter = BackgroundPlotter()
    plotter.add_mesh(deformed_grid, 
                     scalars=colors,rgb=True,
                     scalar_bar_args={"title": "UX"}, 
                     show_edges=True)
    if node_labels:
        plotter.add_point_labels(deformed_grid.points, 
                             labels, font_size=10, point_color='red')
    n_elements = mesh.elements.n_elements        
    plotter.add_text(f"""{file}, {n_points} nodes, {n_elements} elements
scale={scale:.3g}
""", font_size=12)
    return (plotter,scale)
#%% debug cell
qtplot(r_bt, node_labels=True)
#%% st1
# uses cerig which is not compatible with large displacments
if Model.SOLID in models or True:
    do_nlgeom=False
    target_nodes = 110_000
    tolerance = 0.10  # ±10 %
    max_iter = 20
    if target_nodes <= 1000:
        esize=0.12
    elif target_nodes <=30_000:
        esize = 0.02
    elif target_nodes <=120_000:
        esize=0.015
    for i in range(max_iter):
        mapdl.clear()
        mapdl.prep7()
        mapdl.et(1, 187)  # SOLID187 (quadratic tetrahedron)
        mapdl.mp('EX', 1, E)
        mapdl.mp('PRXY', 1, 0.3)
        mapdl.block(0, L, 0, w, 0, h)    
        mapdl.block(0, L,
                    t, w - t,
                    t, h - t)
        mapdl.vsbv(1, 2)
        mapdl.esize(esize)
        mapdl.vmesh('ALL')
        # Get number of elements
        nnodes= int(mapdl.get_value("NODE",0, "COUNT"))
        print(f"Iteration {i+1}: esize={esize:.3g}, nodes={nnodes}")
        # Check if within tolerance
        if abs(nnodes - target_nodes) / target_nodes <= tolerance:
            print("Mesh size target reached. Starting solution")
            break    
        # Update esize based on error
        scale = np.sqrt((target_nodes / nnodes))
        esize /= scale
    else:
        raise Exception("Maximum iterations reached "
                        "without hitting target.")
    mapdl.nsel('S', 'LOC', 'X', 0)
    mapdl.d('ALL', 'ALL', 0)
    cerig_master_node=mapdl.n(x=L, y=w/2, z=h/2)
    mapdl.et(2,'MASS21')
    mapdl.type(2)
    mapdl.tshap('POINT')
    mapdl.r(1)
    mapdl.e(cerig_master_node)
    mapdl.nsel('S', 'LOC', 'X', L)
    mapdl.cm('free_end', 'NODE')
    mapdl.run('CMSEL, S, free_end')
    mapdl.run(f'CERIG, {cerig_master_node}, ALL, UY, UZ')
    mapdl.f(cerig_master_node,'MX',moment)
    mapdl.allsel()
    r_st1=pick_results(mapdl,file='st1')
    r_st1.rfe=r_st1.displacement[cerig_master_node-1][3]*180/np.pi
    print(("Rotation at free end using solid with nlgeom=off (st1)"
           f"{r_st1.rfe:.4g}°"))
    if do_nlgeom:
        r_st1_nl=pick_results(mapdl,True,'st1_nl')
        r_st1_nl.rfe=r_st1_nl.displacement[cerig_master_node-1][3]*180/np.pi
        print(("Rotation at free end using solid with nlgeom=on (st1_nl)"
               f"{r_st1_nl.rfe:.4g}°"))
    if do_plots:
        dpf1(r_st1)
        if do_nlgeom:
            dpf1(r_st1_nl)
        dpf2(r_st1)
        if do_nlgeom:
            dpf2(r_st1_nl)
#%% st2
# uses ce:s which allow shrinking and expansion of cross section
if Model.SOLID in models:
    mapdl.allsel()
    r_st2=pick_results(mapdl,file='st2')
    r_st2_nl=pick_results(mapdl,True,'st2_nl')   
    if do_plots:
        dpf1(r_st2)
        dpf1(r_st2_nl)
        dpf2(r_st2)
        dpf2(r_st2_nl)     
#%% plot results
def get_sorted_node_numbers(result):
    nnum=result.mesh.nnum
    nodes=result.mesh.nodes
    ss=nnum.tolist()
    ss.sort(key=lambda c:nodes[c-1][0])
    return ss

def plot_result(fig,ax,result,index,**kwargs):
    nnum =result.nnum
    disp = result.displacement  # disp shape: (n_nodes, 7)
    # Get coordinates
    coords = result.coords  # shape: (n_nodes, 3)
    # Extract X and ROTX (index 3), convert to degrees
    val = []
    scaler=1
    if rotation_in_degress and index in(3,4,5):
        scaler=(180 / np.pi);
    for i, node in enumerate(nnum):
        x = coords[i][0]
        dval = disp[i][index]
        val.append((x,scaler*dval))
    # Sort by X
    sorted_val = sorted(val, key=lambda pair: pair[0])
    x_vals, vals = zip(*sorted_val)
    # Plot
    ax.plot(x_vals, vals, **kwargs)

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

def add_analytical_rotation(ax,It,Iw):
    if 'r_bt' in globals():
        xv=np.sort(r_bt.coords[:,0])
    else:
        xv=np.linspace(0,L)
    yv=theta(moment,It,Iw,L,xv)
    if rotation_in_degress:
        yv=180/np.pi*yv
    ax.plot(xv,yv,label='analytical')

secdata=bm()
if force:    
    fig_hf, ax_hf = plt.subplots(num='horizontal force',clear=True)
    ax_hf.set_xlabel(r'x-coordinate [m]')
    ax_hf.set_ylabel(r'horizontal displacement [m]')
    if 'r_bhf' in vars():
        plot_result(fig_hf,ax_hf,r_bhf,1,'beam188')
    if 'r_bhf_nl' in vars():
        plot_result(fig_hf,ax_hf,r_bhf_nl,1,'beam188-nlgeom')
    add_analytical_bending(ax_hf,get_sec_property('Izz'))
    ax_hf.legend()
if force:    
    fig_vf, ax_vf = plt.subplots(num='vertical force',clear=True)
    ax_vf.set_xlabel(r'x-coordinate [m]')
    ax_vf.set_ylabel(r'vertical displacement [m]')
    if 'r_bvf' in vars():
        plot_result(fig_vf,ax_vf,r_bvf,2,'beam188')
    if 'r_bvf_nl' in vars():
        plot_result(fig_vf,ax_vf,r_bvf_nl,2,'beam188-nlgeom')
    add_analytical_bending(ax_vf,get_sec_property('Iyy'))
    ax_vf.legend()
if moment:
    fig_t, ax_t = plt.subplots(num='torsion',clear=True)
    ax_t.set_xlabel(r'x-coordinate [m]')
    if rotation_in_degress:
        ax_t.set_ylabel(r'rotation [degrees]')
    else:
        ax_t.set_ylabel(r'rotation [radians]')
    if 'r_bt' in vars():
        plot_result(fig_t,ax_t,r_bt,3
                    ,label='beam188(bt)'
                    ,marker=MarkerStyle((3,0,0))
                    ,markevery=(0.02,0.2)
                    )
    if 'r_bt_nl' in vars():
        plot_result(fig_t,ax_t,r_bt_nl,3
                    ,label='beam188-nlgeom(bt)'
                    ,marker=MarkerStyle((5,0,0))
                    ,markevery=(0.05,0.2)
                    )
    if 'r_bt1' in vars():
        plot_result(fig_t,ax_t,r_bt1,3
                    ,label='beam188(bt1)'
                    ,marker=MarkerStyle((3,1,0))
                    ,markevery=(0.08,0.2)
                    )
    if 'r_bt1_nl' in vars():
        plot_result(fig_t,ax_t,r_bt1_nl,3
                    ,label='beam188-nlgeom(bt1)'
                    ,marker=MarkerStyle((5,1,0))
                    ,markevery=(0.11,0.2)
                    )
    if 'r_bt2' in vars():
        plot_result(fig_t,ax_t,r_bt2,3
                    ,label='beam188(bt2)'
                    ,marker=MarkerStyle((3,2,0))
                    ,markevery=(0.14,0.2)
                    )
    if 'r_bt2_nl' in vars():
        plot_result(fig_t,ax_t,r_bt2_nl,3
                    ,label='beam188-nlgeom(bt2)'
                    ,marker=MarkerStyle((5,2,0))
                    ,markevery=(0.17,0.2)
                    )
    if 'r_st1' in vars() and False:
        plot_solid_result(fig_t,ax_t,r_st1,5,'solid187-cerig(st1)')
    if 'r_st1_nl' in vars() and False:
        plot_solid_result(fig_t,ax_t,r_st1_nl,5,'solid187-cerig(st1)')
    add_analytical_rotation(ax_t,
                           get_sec_property(secdata,'Torsion Constant'),
                           get_sec_property(secdata,'Warping Constant'))
    ax_t.legend()
#%%