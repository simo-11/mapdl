# mapdl
Structural analysis using Ansys Mechanical Parametric Design Language.

# Target
 * Provide resusable reference case for other areas of my interests

# starting grpc server
 * use suitable work directory which is listed in .gitignore
```
C:\Users\simon\github\mapdl\wrk [main ≡]> start-process -FilePath 'C:\Program Files\ANSYS Inc\ANSYS Student\v242\ANSYS\bin\winx64\ANSYS242.exe' -ArgumentList "-grpc" -NoNewWindow
```
Port is shown in output window e.g. Server listening on : 0.0.0.0:50052

# Using spyder
Tools/Preferences/Python interpreter github/mapdl/.venv/Scripts/python.exe
```
In: from ansys.mapdl import core as pymapdl
In: mapdl = pymapdl.Mapdl()
```

# References
 * https://www.ansys.com/blog/what-is-apdl
 * https://mapdl.docs.pyansys.com/
