import numpy as np
import matplotlib.pyplot as plt

data2 = np.loadtxt("data2.dat")


mean = np.mean(data2, axis=0)
std = np.std(data2, axis=0)

plt.subplot(2, 1, 1)
plt.plot(
    data2[:, 0], data2[:, 1], label=f"KE: {mean[1]:.4}" + "$\pm$" + f"{std[1]:.3}", lw=3
)
plt.plot(
    data2[:, 0], data2[:, 2], label=f"PE: {mean[2]:.4}" + "$\pm$" + f"{std[2]:.3}", lw=3
)
plt.plot(
    data2[:, 0],
    data2[:, 3],
    label=f"Energy: {mean[3]:.4}" + "$\pm$" + f"{std[3]:.3}",
    lw=3,
)

plt.xticks(np.arange(0, 17, 1))
plt.grid()
plt.ylabel("Energy", fontsize=18)
plt.title("Energy Profile", fontsize=14)


plt.legend(prop={"size": 20})

plt.subplot(2, 1, 2)
mom2 = np.loadtxt("momt2.dat")

mom_total = np.sqrt(mom2[:, 1] ** 2 + mom2[:, 2] ** 2 + mom2[:, 3] ** 2)

plt.plot(mom2[:, 0], mom2[:, 1], label="$p_x$", lw=2)
plt.plot(mom2[:, 0], mom2[:, 2], label="$p_y$", lw=2)
plt.plot(mom2[:, 0], mom2[:, 3], label="$p_z$", lw=2)
plt.plot(mom2[:, 0], mom_total, label="|p|", lw=2)
plt.title("Momentum conservation ($10^{-15}$ scale)", fontsize=14)
plt.xticks(np.arange(0, 17, 1))
plt.xlabel("No of iterations*200 (time)", fontsize=20)
plt.ylabel("Momentum", fontsize=18)
plt.grid()
plt.suptitle(
    "MD simulation with 2000 particles and thermostat ($r_c=2.5\sigma$)", fontsize=22
)
plt.legend(prop={"size": 15})
plt.show()
