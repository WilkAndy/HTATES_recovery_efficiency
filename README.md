# Recovery efficiency in high-temperature aquifter thermal energy storage systems

The [PorousFlow module](https://mooseframework.inl.gov/modules/porous_flow/index.html) of the finite-element software [MOOSE](https://github.com/idaholab/moose) was used to simulate water flow in a high-temperatue aquifer thermal energy storage system, in order to compute the recovery efficiency.  This repository holds files associated with these simulations and computations.

## Meshes and input files

MOOSE input files are:

- `Simulations/Generic/GeoTES_radial_KT.i`
- `Simulations/Generic/GeoTES_radial_KT_gmsh.i`

Meshes and the scripts used to generate them can be found in `Simulations/Generic/meshes`.

Annual cycles were run using the input file `GeoTES_radial_KT.i`, which generates a quad mesh with radial refinement within 2 times the thermal radius of the well screen.

Daily cycles were run using the input file `GeoTES_radial_KT_gmsh.i`, which uses meshes generated using Gmsh. The meshes comprise a quad mesh within 2 times the thermal radius of the well screen and extending 5 m into the caps, and a triangular mesh elsewhere. This was done in order to achieve the necessary refinement to capture behaviour within the relatively small thermal radius of the daily cycles, while avoiding the need for such refinement outside this area.

## Water properties

During the course of this study, it was found that the water properties (chiefly density and viscosity) impact the recovery efficiently greatly, so accurate estimation of these is necessary for accurate estimation of the recovery efficiency.

- `Simulations/Generic/water97_tabulated_modified.csv` are properties using the IAPWS-97 equation of state
- `Simulations/Generic/waterBuscheck_tabulated.csv` are properties using Buscheck's approximation
- `Simulations/Generic/waterSchout_tabulated_thermExp2E-4.csv` are properties using Schout's approximation

## Results

Results are in in `Simulations/Generic/REpaper/results`.

- `Simulations/Generic/REpaper/results/RE/all_RE_zero_dispersion.csv` holds a summary of the results assuming zero dispersion
- `Simulations/Generic/REpaper/results/RE/all_RE.csv` holds a summary of all results
- `Simulations/Generic/REpaper/results/csv/` is a directory that holds all the raw MOOSE output, organised by cycle length, depth and dispersion




