# mapdl
Structural analysis using Ansys Mechanical Parametric Design Language through python interface.

# Target
 * Provide resusable reference cases for other areas of my interests

# Update of python packages
```
C:\Users\simon\github\mapdl [main ≡]> .venv\Scripts\activate
(mapdl) C:\Users\simon\github\mapdl [main ≡]> uv pip install ansys-mapdl-core[graphics] -U
(mapdl) C:\Users\simon\github\mapdl [main ≡]> uv pip install ansys-dpf-core -U
(mapdl) C:\Users\simon\github\mapdl [main ≡]> uv pip install ansys-dpf-post[plotting] -U
```
# starting grpc server
grpc version can be started within python. In that case ansys files like file.err and other will be %TEMP%\ansys_...
 * use suitable work directory which is listed in .gitignore e.g. C:\Users\simon\github\mapdl\wrk
```
 start-process -FilePath 'C:\Program Files\ANSYS Inc\ANSYS Student\v242\ANSYS\bin\winx64\ANSYS242.exe' -ArgumentList "-grpc" -NoNewWindow
```
Port is shown in output window e.g. Server listening on : 0.0.0.0:50052

# Using spyder
Tools/Preferences/Python interpreter github/mapdl/.venv/Scripts/python.exe
```
In: from ansys.mapdl import core as pymapdl
In: mapdl = pymapdl.Mapdl()
or
In: mapdl=pymapdl.launch_mapdl()
```

# References
 * https://www.ansys.com/blog/what-is-apdl
 * https://mapdl.docs.pyansys.com/
 * https://math.docs.pyansys.com/
 * https://dpf.docs.pyansys.com/
 * https://post.docs.pyansys.com/version/stable/
 * https://ansyshelp.ansys.com/public/account/secured?returnurl=////Views/Secured/prod_page.html?pn=Mechanical%20APDL&prodver=24.2&lang=en
 * https://docs.pyvista.org/
