# Comparisons with Auburn experimental data

Instructions:
- the run.sh script contains all the parameters and explanations of where in Buscheck's PhD thesis they may be found.  Simply run it.  It copies ../../GeoTES_radial_KT.i and creates the water properties using generate_auburn_water.py, then runs MOOSE
- generate_auburn_water.py generates the tabulated water properties based on Buscheck's PhD thesis
- plot_temperature.py plots the results, which are assumed to be in gold/auburn_minmod.csv


Buscheck's PhD thesis~\cite{buscheck1984} carefully describes the second set of ATES experiments at Auburn University~\cite{molz1978,molz1979}.  Hot water was injected into a shallow confined aquifer, stored and then produced from the aquifer, while a large number of temperature observations were performed.  Buscheck presents fully-parameterised conceptual and numerical models.  Using Buscheck's parameters --- the geometric parameters associated with the borehole, aquifer and cap rocks, the physical parameters (porosity, permeability, etc) of the aquifer and cap rocks, the {\em in-situ} temperature, the injection temperature, the average injection and production rates, and the temperature-dependent fluid properties --- in our model, the simulated temperature of the produced hot water may be compared with the observed value.

Figure [auburn_stage1.png] shows that the MOOSE model agrees well with the observations at Auburn University.  The descrepancy at the beginning of the production period is due to the experimental injection temperature and rate during the final days of injection being higher than average.  In contrast, the MOOSE model is simplified in that it does not precisely account for the hourly fluctuations in injection rate and temperature.  After a short time, the impact of those fluctuations disappears and MOOSE matches the experimental observations, demonstrating that the MOOSE model can accurately simulate ATES systems.

\caption{The observed and simulated first-cycle production temperatures as a function of time (observed values sourced from~\cite[Figure III.12]{buscheck1984}).  The descrepancy at the beginning of the production period is due to the experimental injection temperature and rate during the final days of injection being higher than average.}

