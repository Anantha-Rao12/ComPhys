import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt("data1.dat")

mean = np.mean(data1, axis=0)
std = np.std(data1, axis=0)

plt.plot(
    data1[:, 0], data1[:, 1], label=f"KE: {mean[1]:.4}" + "$\pm$" + f"{std[1]:.3}", lw=3
)
plt.plot(
    data1[:, 0], data1[:, 2], label=f"PE: {mean[2]:.4}" + "$\pm$" + f"{std[2]:.3}", lw=3
)
plt.plot(
    data1[:, 0],
    data1[:, 3],
    label=f"Energy: {mean[3]:.4}" + "$\pm$" + f"{std[3]:.3}",
    lw=3,
)

plt.xticks(np.arange(0, 17, 1))
plt.grid()
plt.xlabel("No of iterations (time)", fontsize=18)
plt.ylabel("Energy", fontsize=18)
plt.title("MD simulation with 1000 particles ($r_c=2.5\sigma$)", fontsize=22)

plt.legend(prop={"size": 20})
plt.show()
