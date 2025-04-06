# Initial installation

Based on https://mapdl.docs.pyansys.com/version/stable/getting_started/index.html#ref-getting-started
 * Virtual env using vu from [astral-sh|https://docs.astral.sh/uv/pip/environments/]
```
C:\Users\simon\github\mapdl [main ≡]> uv venv --python 3.12.9
C:\Users\simon\github\mapdl [main ≡]> .venv\Scripts\activate
```
 * install ansys-mapdl-core
```
(mapdl) C:\Users\simon\github\mapdl [main ≡]> uv pip install ansys-mapdl-core
```
 * verify
```
>>> from ansys.mapdl import core as pymapdl
>>> mapdl=pymapdl.launch_mapdl()
>>> print(mapdl)
Mapdl
-----
PyMAPDL Version:     0.69.3
Interface:           grpc
Product:             Ansys Mechanical Enterprise Academic Student
MAPDL Version:       24.2
Running on:          localhost
                     (127.0.0.1)
>>> mapdl.exit()
```
 ** Started instance is shown in task manager like this  ![image](https://github.com/user-attachments/assets/429461d3-5bb1-4c9a-8858-f9d1ecc19cad)

# Update python
 
# Update only python modules
```
C:\Users\simon\github\mapdl [main ≡]> .venv\Scripts\activate
(mapdl) C:\Users\simon\github\mapdl [main ≡]> uv pip install --upgrade ansys-mapdl-core
(mapdl) C:\Users\simon\github\mapdl [main ≡]> deactivate
```
