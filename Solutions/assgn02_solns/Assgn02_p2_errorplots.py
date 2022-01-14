import numpy as np
import matplotlib.pyplot as plt

trap = np.loadtxt("p2_trap_e3.dat")
mc_sample = np.loadtxt("p2_MC_sampling.dat")
mc_hitmiss = np.loadtxt("exp_accep_rej.dat")

data = [trap, mc_sample, mc_hitmiss]
labels = [
    "Trapezium method",
    "Monte Carlo Acceptance rejection method",
    "Monte Carlo importance sampling",
]
# Plotting data

markers = ["o", "p", "D"]
for idx, dataset in enumerate(data):

    plt.plot(
        np.log10(dataset[:, 0]),
        dataset[:, 1],
        marker=markers[idx],
        alpha=0.5,
        markersize=10,
        markeredgewidth=5,
        ls="None",
        label=labels[idx],
    )
    plt.errorbar(
        x=np.log10(dataset[:, 0]),
        y=dataset[:, 1],
        ls="None",
        yerr=dataset[:, 2],
        capsize=5,
        color="k",
        elinewidth=2,
    )

plt.axhline(y=np.exp(3) - 1, color="k", ls="--", lw=3, label="Actual Value", alpha=0.5)
plt.xlabel("$Log_{10}$(x)", fontsize=18)
plt.ylabel("Value of the integral", fontsize=18)
plt.title(
    "Comparison of errors from different integration method to solve $\int_0^3 e^x dx$",
    fontsize=20,
)
plt.yticks(np.arange(-20, 50, 10), np.arange(-20, 50, 10), fontsize=14)
plt.xticks(np.arange(0, 10), np.arange(10), fontsize=14)
plt.legend(prop={"size": 15})
plt.show()
