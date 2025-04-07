# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 10:48:02 2025

@author: simo nikula

Create, solve and report results using ansys models for cantilever U-section
"""
# %% settings and commons initializations
from ansys.mapdl import core as pymapdl
h=0.1
w=0.05
t=0.004
E=210E9
nu=0.3
sharp_corners=True
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
mapdl=None
def initMapdl():
    global mapdl
    if mapdl is None:
        mapdl=pymapdl.Mapdl()
initMapdl()        
# %% Beam model
mapdl.clear()
mapdl.prep7()
mapdl.et(1,"BEAM188")
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
mapdl.dk(kpoi=1,lab="ALL")
mapdl.nplot(vtk=True, nnum=True, cpos="xy",
            show_bounds=True, point_size=10)
# %% beam vertical force
# %% beam horizontal force
# %% beam torsion
if moment is None:
    raise Exception(no_moment_defined)

