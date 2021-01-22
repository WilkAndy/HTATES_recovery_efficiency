# HTATES_recovery_efficiency
Files associated with computing recovery efficiency in HT-ATES systems

Meshes and input files:

Annual cycles were run using the input file GeoTES_radial_KT.i, which generates a quad mesh with radial refinement within 2 times the thermal radius of the well screen.

Daily cycles were run using the input file GeoTES_radial_KT_gmsh.i, which uses meshes generated using Gmsh. The meshes comprise a quad mesh within 2 times the thermal radius of the well screen and extending 5 m into the caps, and a triangular mesh elsewhere. This was doen in order to achieve the necessary refinement to capture behaviour within the relatively small thermal radius of the daily cycles, while avoiding the need for such refinement outside this area. The meshes and the scripts used to generate them can be found in Simulations/Generic/meshes.

