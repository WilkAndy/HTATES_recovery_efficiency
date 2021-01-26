#!/usr/bin/python3

# Generates density, viscosity, enthalpy and internal_energy based on the functions presented in Buscheck "The hydrothermal analysis of aquifer themal energy storage" (1984) PhD thesis, University of California Berkeley

import sys

def density(tK):
    # Density, given Temperature in Kelvin.
    # page 35 of Buscheck
    tC = tK - 273.15
    if (tC >= 25):
        return 996.9 * (1 + (tC - 25.0) * (-3.17E-4 + (tC - 25.0) * (-2.56E-6)))
    return 996.9 * (1 - 1.87E-4 * (tC - 25.0))

def viscosity(tK):
    # Viscosity, given Temperature in Kelvin
    # page 35 of Buscheck
    tC = tK - 273.15
    if (tC <= 20.0):
        return 0.1005E-2
    elif (tC <= 50.0):
        return 0.545E-3 + (0.1005E-2 - 0.545E-3) * (tC - 50.0) / (20.0 - 50.0)
    elif (tC <= 100.0):
        return 0.280E-3 + (0.545E-3 - 0.280E-3) * (tC - 100.0) / (50.0 - 100.0)
    elif (tC <= 150.0):
        return 0.182E-3 + (0.280E-3 - 0.182E-3) * (tC - 150.0) / (100.0 - 150.0)
    return 0.182E-3

def internal_energy(tK):
    # Internal energy, given Temperature in Kelvin
    # page 35 of Buscheck
    # Gao et al assume a constant heat capacity (cp), and internal_energy = cp * TC, see Eqn (2)
    tC = tK - 273.15
    if (tC <= 20.0):
        cp = 0.4182E4
    elif (tC <= 75.0):
        cp = 0.3894E4 + (0.4182E4 - 0.3894E4) * (tC - 75.0) / (20.0 - 75.0)
    elif (tC <= 125.0):
        cp = 0.3652E4 + (0.3894E4 - 0.3652E4) * (tC - 125.0) / (75.0 - 125.0)
    elif (tC <= 200.0):
        cp = 0.3341E4 + (0.3652E4 - 0.3341E4) * (tC - 200.0) / (125.0 - 200.0)
    else:
        cp = 0.3341E4
    return cp * tC
    
def enthalpy(tK):
    # Enthalpy, given Temperature in Kelvin
    # page 35 of Buscheck
    return internal_energy(tK)

fn = "waterAuburn_tabulated.csv"
sys.stdout.write("Writing " + fn + "\n")
f = open(fn, "w")
f.write("pressure, temperature, density, viscosity, enthalpy, internal_energy\n")
for p in [0, 5E7]:
    for tK in range(275, 551):
        f.write(str(p) + ", " + str(tK) + ", " + str(density(tK)) + ", " + str(viscosity(tK)) + ", " + str(enthalpy(tK)) + ", " + str(internal_energy(tK)) + "\n")
f.close()
sys.stdout.write("Done\n")
sys.exit(0)
