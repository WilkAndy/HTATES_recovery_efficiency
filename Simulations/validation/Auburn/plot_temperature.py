import matplotlib.pyplot as plt
import sys

f = open("gold/auburn_minmod.csv", "r")
data = [line.strip().split(",") for line in f.readlines()[1:]]
f.close()
hours = [float(d[0]) * 24 for d in data]
temp = [float(d[6]) for d in data]

obs_hours = [3113, 3140, 3281, 3449, 3617, 3785, 3953, 4090]
obs_temp = [56.0,  55.0, 50.0, 46.0, 42.0, 38.3, 35.0, 32.7]

plt.figure()
plt.plot(hours, temp, 'k-', label="Simulated")
plt.plot(obs_hours, obs_temp, "ks", label="Observed")
plt.legend()
plt.gca().set_xlim([3113, 4100])
plt.gca().set_ylim([32, 57])
plt.xlabel("Time (hours)")
plt.ylabel("Temperature (degC)")
plt.title("Temperature of produced water in the Auburn experiment")
plt.grid()
plt.savefig("auburn_stage1.png", bbox_inches="tight")
plt.show()
sys.exit(0)

