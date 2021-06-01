#!/bin/bash -l

# Parameterisation from the Buscheck 1984 PhD thesis
# aq_thickness = 21m (page 28)
# cap_thickness = 9m (page 28)
# screen_top = 10.5 (page 28 - "screened in the upper half", but picture on page 41 looks like about screen_top=9 - not sure which is correct)
# screen_bottom = 0 (page 28 - "screened in the upper half", but picture on page 41 looks like about screen_bottom=-1 - not sure which is correct)
# depth = 20m (page 28 says depth of aquifer top is 40m, but picture on page 36 shows 20m - presumably the water table is 20m below the topography)
# insitu_temp = 20degC (page 28)
# inject_temp = 55.2degC (page 28.  However, note that inject_heat1/function is chosen below to reflect the varying temperature of the injection and the varying injection rate.  This is only an approximate function, estimated from the graph FigureIII.4 on page 40, based on the true injection rate and injection temperature.  The average value of injection temperature = 55.2degC.)
# fluid density = 996.9 * (1 - 3.17E-4 * (T - 25) - 2.56E-6 * (T - 25)^2), where T is temperature in degC (page 35)
# At T = 55.2degC, fluid density = 985kg.m^-3
# injection volume  = 55000m^3 (page 28) , so inject_fluid_mass = 5.4177E7kg
# inject_time = 79.166667day (=1900 hours, page 28)
# store_time = 50.5416667day (1213 hours, page 28)
# produce_time = 41.125day (=987hour, page 28)
# produce_fluid_mass = 5.56076E7kg (=15.65kg/s * 987hr * 3600, page 28)
# aq_porosity = 0.25 (page 34)
# aq_hor_perm = 5.3E-11m^2 (page 32)
# aq_ver_perm = 5.3E-12m^2 (page 33)
# aq_density = 2.6E3kg/m^3 (page 32)
# aq_specific_heat_cap = 696J/kg/K (page 32)
# aq_hor_thermal_cond = 2.29W/m/K (page 32)
# aq_ver_thermal_cond = 2.29W/m/K (page 32)
# cap_porosity = 0.15 (page 34)
# cap_hor_perm = 5.3E-16m^2 (page 33)
# cap_ver_perm = 5.3E-17m^2 (page 33)
# cap_Density = 2.6E3kg/m^3 (page 32)
# cap_specific_heat_cap = 696J/kg/K (page 32)
# cap_hor_thermal_cond = 2.56W/m/K (page 32)
# cap_ver_thermal_cond = 2.56W/m/K (page 32)
echo "Defining directories..."
export COMBINED_DIR=/Users/wil04q/projects/moose/modules/combined

cp ../../GeoTES_radial_KT.i .
python generate_auburn_water.py

echo "Running..."
 mpirun -np 2 ${COMBINED_DIR}/combined-opt -i GeoTES_radial_KT.i --show-input \
	 flux_limiter=minmod  \
	 num_cycles=1  \
	 depth=20  \
	 aq_thickness=21 \
	 screen_top=10.5 \
	 screen_bottom=0 \
	 inject_time=79.1666667  \
	 store_time=50.5416667  \
	 produce_time=41.125  \
	 rest_time=0  \
	 cap_thickness=9  \
	 insitu_temp=20  \
	 geothermal_gradient=0  \
	 aq_porosity=0.25  \
	 aq_hor_perm=5.3E-11 \
	 aq_ver_perm=5.3E-12 \
	 aq_density=2600.0  \
	 aq_specific_heat_cap=696  \
	 aq_hor_thermal_cond=2.29  \
	 aq_ver_thermal_cond=2.29  \
	 cap_porosity=0.15  \
	 cap_hor_perm=5.3E-16  \
	 cap_ver_perm=5.3E-17  \
	 cap_density=2600.0  \
	 cap_specific_heat_cap=696  \
	 cap_hor_thermal_cond=2.56  \
	 cap_ver_thermal_cond=2.56  \
	inject_fluid_mass=5.4177E7 \
	produce_fluid_mass=5.56076E7 \
	inject_temp=55.2 \
	Materials/porosity_caps/biot_coefficient=0.6 \
	Materials/porosity_aq/biot_coefficient=0.6 \
	Materials/porosity_caps/solid_bulk=7.21E7 \
	Materials/porosity_aq/solid_bulk=7.21E7 \
	Outputs/inactive='temperature_on_line' \
	Postprocessors/inactive='' Functions/inactive=''\
	filename='auburn' \
	Modules/FluidProperties/tabulated_water/fluid_property_file=waterAuburn_tabulated.csv \
	Outputs/checkpoint=false \
	Executioner/nl_max_its=30 \
	Executioner/dtmax=1.0 \
	Executioner/TimeStepper/optimal_iterations=5 \
	BCs/inject_heat1/function='if(t<=14,56,if(t<=24,57,if(t<=56,52.85,if(t<=70,56,if(t<=80,59,${inject_temp})))))'
echo "Finished"
