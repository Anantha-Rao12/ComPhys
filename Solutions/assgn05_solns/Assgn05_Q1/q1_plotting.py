import numpy as np
import matplotlib.pyplot as plt

euler = np.loadtxt("q1_eulermethod.dat")
modified_euler = np.loadtxt("q1_modified_euler_method.dat")
improved_euler = np.loadtxt("q1_improved_euler_method.dat")
rk4 = np.loadtxt("q1_rk4.dat")

dataset = [euler, modified_euler, improved_euler, rk4]
legends = [
    "Euler, dx=0.001",
    "Mod. Euler, dx=0.001",
    "Imp. Euler, dx=0.001",
    "RK4, dx=0.01",
]
marker = ["p", "P", "o", "s"]

plt.xlabel("X (=h*niterations)", fontsize=20)
plt.ylabel("Y =tan(x)", fontsize=20)
for id, data in enumerate(dataset):
    plt.plot(
        data[:, 0], data[:, 1], lw=1, label=legends[id], ls="--", marker=marker[id]
    )


plt.plot(euler[:, 0], np.tan(euler[:, 0]), label="Tan(x) (theoretical)")
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.legend(prop={"size": 15})
plt.show()
