# Recovery efficiency in high-temperature aquifter thermal energy storage systems

The PorousFlow module (https://mooseframework.inl.gov/modules/porous_flow/index.html) of the finite-element software MOOSE (https://github.com/idaholab/moose) was used to simulate water flow in a high-temperatue aquifer thermal energy storage system, in order to compute the recovery efficiency.  This repository holds files associated with these simulations and computations.

For full details refer to Sheldon et al. (2021).

## Meshes and input files

MOOSE input files are:

- Simulations/GeoTES_radial_KT.i
- Simulations/GeoTES_radial_KT_gmsh.i

Note that these input files may not work with future versions of MOOSE due to syntax changes. An updated version can be found in the MOOSE GitHub repository, located in: https://github.com/idaholab/moose/tree/master/modules/porous_flow/examples/ates

Annual cycles were run using the input file GeoTES_radial_KT.i, which generates a quad mesh with radial refinement within 2 times the thermal radius of the well screen.

Daily cycles were run using the input file GeoTES_radial_KT_gmsh.i, which uses meshes generated using Gmsh. Meshes and the scripts used to generate them can be found in Simulations/meshes. The meshes comprise a quad mesh within 2 times the thermal radius of the well screen and extending 5 m into the caps, and a triangular mesh elsewhere. This was done in order to achieve the necessary refinement to capture behaviour within the relatively small thermal radius of the daily cycles, while avoiding the need for such refinement outside this area.

These input files were modified from the command line to reflect the aquifer properties and operating parameters for each scenario.

## Water properties

During the course of this study, it was found that the water properties (chiefly density and viscosity) impact the recovery efficiently greatly, so accurate estimation of these is necessary for accurate estimation of the recovery efficiency.

- Simulations/water97_tabulated_modified.csv are properties obtained using the IAPWS-97 equation of state
- Simulations/waterBuscheck_tabulated.csv are properties obtained from the approximation used by Buscheck (1984)
- Simulations/waterSchout_tabulated_thermExp2E-4.csv are properties obtained from the approximation used by Schout et al. (2014)

## Results

Results are in Simulations/results.

- Simulations/results/RE/all_RE_zero_dispersion.csv holds a summary of the results assuming zero dispersion
- Simulations/results/RE/all_RE.csv holds a summary of all results
- Simulations/results/csv/ is a directory that holds all the raw MOOSE output, organised by cycle length, depth and dispersion

## Correlations

The architecture, training and testing of a simple convolutional neural network for predicting recovery efficiency from aquifer and operating parameters is defined through the python script CNN/plot_nn1.py.

## Validation

The model was validated by comparison with 3 previous studies:

1. Reproducing the results of a field experiment conducted by Auburn university (files in Simulations/validation/Auburn);
2. Reproducing some numerical models by Buscheck (1984) (files in Simulations/validation/Buscheck);
3. Reproducing some numerical models by Schout et al. (2014) (files in Simulations/validation/Schout). 

## References

Buscheck, T.A.A., 1984. The hydrothermal analysis of aquifer themal energy storage. PhD Thesis, University of California, Berkeley.

Schout, G., Drijver, B., Gutierrez-Neri, M., Schotting, R., 2014. Analysis of recovery efficiency in high-temperature aquifer thermal energy storage: a Rayleigh-based method. Hydrogeol. J. 22, 281–291. https://doi.org/10.1007/s10040-013-1050-8

Sheldon, H.A., Wilkins, A. and Green, C.P., 2021 (in review). Recovery efficiency in high-temperature aquifer thermal energy storage systems. Geothermics.
