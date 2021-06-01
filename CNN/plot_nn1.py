# Using all_RE_zero_dispersion.csv, constructs two 'best fit' CNNs: one for each injection time.  Plots the results
from numpy.random import seed
seed(1)
from tensorflow import set_random_seed
set_random_seed(2)

import os
import sys
from bisect import bisect_left
import math
import matplotlib.pyplot as plt

import numpy as np
from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import r2_score

def waterProps(p, TC, data, p_vals, T_vals):
    """
    returns water properties at (p, TC) using "data", which is assumed to be tabulated in the form of water97_tabulated_modified.csv
    @param p pressure in Pascals
    @param TC temperature in degC
    @param data (p, TK, density, viscosity, enthalpy, internal_energy) read from water97_tabulated_modified.csv
    @param p_vals sorted pressure values in water97_tabulated_modified.csv
    @param T_vals sorted temperature (Kelvin) values in water97_tabulated_modified.csv
    """
    p_low_index = bisect_left(p_vals, p)
    if (p_low_index == 0):
        sys.stderr.write("Pressure " + str(p) + " too low\n")
        sys.exit(1)
    p_low = p_vals[p_low_index - 1]
    p_high = p_vals[p_low_index]

    TK = TC + 273.15
    t_low_index = bisect_left(T_vals, TK)
    if (t_low_index == 0):
        sys.stderr.write("Temperature " + str(TC) + " too low\n")
        sys.exit(1)
    T_low = T_vals[t_low_index - 1]
    T_high = T_vals[t_low_index]

    low_low = [line for line in data if line[0] == p_low and line[1] == T_low][0]
    low_high = [line for line in data if line[0] == p_low and line[1] == T_high][0]
    high_low = [line for line in data if line[0] == p_high and line[1] == T_low][0]
    high_high = [line for line in data if line[0] == p_high and line[1] == T_high][0]

    interpolated = [(((low_low[i] * (p_high - p) + high_low[i] * (p - p_low)) / (p_high - p_low)) * (T_high - TK) + ((low_high[i] * (p_high - p) + high_high[i] * (p - p_low)) / (p_high - p_low)) * (TK - T_low)) / (T_high - T_low) for i in range(6)]
    dint_dp = [(((low_low[i] * (-1) + high_low[i] * (1)) / (p_high - p_low)) * (T_high - TK) + ((low_high[i] * (-1) + high_high[i] * (1)) / (p_high - p_low)) * (TK - T_low)) / (T_high - T_low) for i in range(6)]
    dint_dT = [(((low_low[i] * (p_high - p) + high_low[i] * (p - p_low)) / (p_high - p_low)) * (-1) + ((low_high[i] * (p_high - p) + high_high[i] * (p - p_low)) / (p_high - p_low)) * (1)) / (T_high - T_low) for i in range(6)]
    density = interpolated[2]
    viscosity = interpolated[3]
    enthalpy = interpolated[4]
    internal_energy = interpolated[5]
    bulk_mod = density / dint_dp[2]
    thermal_exp = - dint_dT[2] / density
    vol_specific_heat = dint_dT[5] * density

    return (density, viscosity, enthalpy, internal_energy, bulk_mod, thermal_exp, vol_specific_heat)


# Read definitions and results from all simulations (ignore first column, which is Scenario name)
f = open("../Simulations/results/RE/all_RE_zero_dispersion.csv", 'r', encoding='utf-8-sig')
header = [h.strip() for h in f.readline().strip().split(",")[1:]]
data = [list(map(float, line.strip().split(",")[1:])) for line in f.readlines()]
f.close()

# figure out the column numbers for each property or results
col = {} # col["name"] is the column for "name" in the file
for i in range(len(header)):
    col[header[i]] = i

desired_data = ["perm_x_aq", "perm_y_aq", "perm_x_cap", "perm_y_cap", "depth", "aq_thickness", "injected_fluid_mass", "injected_temp", "inject_time", "store_time", "produce_time", "rest_time", "T_ambient", "P_ambient", "enthalpy_i", "enthalpy_a", "RE1", "RE2", "RE3", "RE4", "RE5"]

# check data is as desired
for nm in desired_data:
    if not (nm in col):
        sys.stderr.write("The column " + nm + " does not appear in " + fn + "\n")
        sys.exit(1)

# extract the thicknesses and injection times simulated
thicknesses = sorted(list(set([line[col["aq_thickness"]] for line in data])))
inject_times = sorted(list(set([line[col["inject_time"]] for line in data])))
sys.stdout.write("Simulations performed for thicknesses " + " ".join(map(str, thicknesses)) + "\n")
sys.stdout.write("Simulations performed for injection times " + " ".join(map(str, inject_times)) + "\n")


# Read the water97 properties
f = open("../Simulations/water97_tabulated_modified.csv", "r")
water97 = [list(map(float, line.strip().split(","))) for line in f.readlines()[1:]]
f.close()
p_vals = [x[0] for x in water97]
seen = set();
p_vals = [x for x in p_vals if x not in seen and not seen.add(x)]
T_vals = [x[1] for x in water97]
seen = set();
T_vals = [x for x in T_vals if x not in seen and not seen.add(x)]


# Compute the modified Rayleigh number
sys.stdout.write("Computing modified Rayleigh\n")
for i in range(len(data)):
    line = data[i]
    T_i = line[col["injected_temp"]]
    T_a = line[col["T_ambient"]]
    H_a = line[col["aq_thickness"]]
    k_ah = line[col["perm_x_aq"]]
    k_av = line[col["perm_y_aq"]]
    d_a = line[col["depth"]]
    p_a = line[col["P_ambient"]]
    (density_amb, viscosity_amb, enthalpy_amb, internal_energy_amb, bulk_mod_amb, thermal_exp_amb, vol_specific_heat_amb) = waterProps(p_a, T_a, water97, p_vals, T_vals)
    (density_i, viscosity_i, enthalpy_i, internal_energy_i, bulk_mod_i, thermal_exp_i, vol_specific_heat_i) = waterProps(p_a, T_i, water97, p_vals, T_vals)
    porosity = 0.25
    heat_cap_a = 800 * 2650
    C_a = porosity * vol_specific_heat_i + (1.0 - porosity) * heat_cap_a
    R_th = math.sqrt(line[col["injected_fluid_mass"]] / density_i * vol_specific_heat_i / (C_a * math.pi * H_a))
    lambda_a = 3 # no dispersion
    Ra_star = thermal_exp_i * (0.5 * (density_amb + density_i)) * 9.81 * pow(H_a, 2) * C_a * math.sqrt(k_av * k_ah) * (T_i - T_a) / (0.5 * (viscosity_amb + viscosity_i)) / lambda_a / R_th
    data[i].append(Ra_star)
len_col = len(col)
col["Ra_star"] = len_col

def cnn_model():
	# create model.  Note: input_dim and other numbers depend on the number of predictors (columns in x)
	model = Sequential()
        # Just a single, fully connected, hidden layer with XXXX neurons
	model.add(Dense(12, input_dim=8, kernel_initializer='normal', activation='relu'))
	model.add(Dense(1, kernel_initializer='normal'))
	# Compile model
	model.compile(loss='mean_squared_error', optimizer='adam')
	return model

fig = plt.figure(0)
styles = ['k', 'r', 'g', 'y', 'c']
for inj_ind in range(len(inject_times)):
    inject_time = inject_times[inj_ind]
    x = [[math.log10(line[col["perm_x_aq"]]), math.log10(line[col["perm_y_aq"]]), line[col["depth"]], line[col["aq_thickness"]], line[col["injected_temp"]], line[col["T_ambient"]], line[col["enthalpy_a"]], line[col["Ra_star"]]] for line in data if (line[col["inject_time"]] == inject_time)]
    #x = [[math.log10(line[col["perm_x_aq"]]), math.log10(line[col["perm_y_aq"]]), line[col["depth"]], line[col["aq_thickness"]], line[col["injected_temp"]], line[col["Ra_star"]]] for line in data if (line[col["inject_time"]] == inject_time)]
    y = np.array([line[col["RE5"]] for line in data if (line[col["inject_time"]] == inject_time)])
    # normalise data
    means = [sum([v[i] for v in x]) / len(x) for i in range(len(x[0]))]
    sds = [np.sqrt(sum([pow(v[i] - means[i], 2) / len(x) for v in x])) for i in range(len(x[0]))]
    x = np.array([np.array([(x[j][i] - means[i]) / sds[i] for i in range(len(x[0]))]) for j in range(len(x))])

    # CNN
    estimator = KerasRegressor(build_fn=cnn_model, epochs=50, batch_size=5, verbose=1) # define the parameters of the CNN
    n_splits = 10 # TODO
    kfold = KFold(n_splits = n_splits, shuffle = True) # define how the (x,y) data will be split into training and test data.
    predicted_re = []
    actual_re = []
    all_inject_time = []
    r2 = []
    fld = -1
    for train_index, test_index in kfold.split(x):
        fld += 1
        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        estimator.fit(x_train, y_train)
        predicted_re += estimator.predict(x_test).tolist()
        actual_re += y_test.tolist()
        all_inject_time += [v[-1] for v in x[test_index].tolist()]
        r2.append(r2_score(y_test, estimator.predict(x_test)))
        sys.stdout.write("  Fold " + str(fld) + " of " + str(n_splits) + ".  R-squared on RE predictions = " + str(r2[-1]) + "\n")
    mean_actual = sum(actual_re) / len(actual_re)
    ss_actual = sum([pow(a - mean_actual, 2) for a in actual_re])
    ss_error = sum([pow(actual_re[i] - predicted_re[i], 2) for i in range(len(predicted_re))])
    av_r2 = str(round(1.0 - ss_error / ss_actual, 2))
    #av_r2 = str(round(sum(r2) / len(r2), 2))
    sys.stdout.write("  Average R-squared on RE predictions = " + av_r2 + "\n")

    # plot the results
    label = "annual cycle"
    if inject_times[inj_ind] == 0.3333:
        label = "daily cycle"
    plt.plot(actual_re, predicted_re, styles[inj_ind] + '.', label = label + "  $R^{2}$ = " + av_r2)
    sys.stdout.write("  R-squared on RE predictions for " + label + " = " + av_r2 + "\n")
plt.plot([0, 1], [0, 1], 'k-', label = "Perfect agreement")
plt.grid()
plt.xlabel("Actual $R$")
plt.ylabel("Predicted $R$")
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.legend()
plt.title("Actual recovery efficiency vs predicted using CNN")
plt.savefig("nn1_correlation.png")
plt.show()

sys.exit(0)
