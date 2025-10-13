# IMAP

IMAS mapping framework written in Python.
Mapping configuration files are written in Jinja-templated YAML files.

## Dependencies

- IMAS-Python
- MDSplus Python interface (or mdsthin)
- ruamel.yaml (YAML 1.2)
- Jinja2
- FastAPI, Uvicorn (Optional, required to run ASGI server)

## List of IDSs available (ongoing)

- charge_exchange
- coils_non_axisymmetric
- ec_launchers
- ece
- equilibrium (use equilibrium/postprocess to extend calculations)
- interferometer
- magnetics
- nbi
- pf_active
- tf
- thomson_scattering
- wall

## Usage

### Mappings

- Mapping equilibrium and wall IDSs for shot #36996
  ```
  python mapping.py -s 36996 equilibrium wall
  ls 36996/ # prints equilibrium.h5  master.h5  wall.h5
  ```
- Mapping ec_launchers and nbi IDSs for shot #36996, during time interval [0.0, 10.0)
  ```
  python mapping.py ec_launchers nbi -s 36996 --tbegin 0.0 --tend 10.0
  ```
- Mapping equilibrium (with or without post-processing) with different trees
  ```
  python mapping.py -s 36996 "equilibrium?postprocess=true&tree=efitrt1" wall -p ./36996_efitrt1
  python mapping.py -s 36996 "equilibrium?postprocess=false" wall -p ./36996_efit01
  python mapping.py -s 36996 "equilibrium?tree=efit02" wall -p ./36996_efit02
  ```

### Plotting and checking

Please check more plotting scripts in `./plotting` directory

- Plotting 2d/3d machine
  ```
  python plotting/plot_2d.py ./36996
  python plotting/plot_3d.py ./36996
  ```
- Compare equilibrium mappings
  ```
  python plotting/equilibrium.py ./36996_efitrt1 ./36996_efit01 ./36996_efit02 3.0
  ```
- Plotting nbi or ec_launchers mapping
  ```
  python plotting/nbi.py ./36996
  python plotting/ec_launchers.py ./36996
  ```
- Checking mapped data
  ```python
  import imas
  db = imas.DBEntry("imas:hdf5?path=./36996", 'r') # data directory
  wall = db.get("wall") # mapped IDS name
  wall.description_2d[0].limiter.unit[0].outline.r # prints array([1.265, ...])
  ```

### Server/client

- Running server
  ```
  python server.py
  ```
  
- Fetching data on the client side
  ```python
  import requests
  response = requests.post("http://localhost:8000/equilibrium", json={"shot": 36996})

  import imas
  factory = imas.IDSFactory("3.39.0")
  equilibrium = factory.equilibrium()
  equilibrium.deserialize(response.content)
  print(equilibrium.time) # prints array of floats if successful
  ```

## Environment settings

### TRANSPNODE
```
module load mdsplus
module load gcc
module load conda3/python3.11.5-al-dd3.39.0
. /XR_TRANSP/app/imas-al/bin/al_env.sh
```
