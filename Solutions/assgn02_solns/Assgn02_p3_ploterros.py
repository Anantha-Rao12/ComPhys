import numpy as np
import matplotlib.pyplot as plt

bruteforce = np.loadtxt("multi_brute_mc.dat")
isample = np.loadtxt("multi_impsampling_mc.dat")

dataset = [bruteforce, isample]
markers = ["o", "P"]
labels = ["Monte-Carlo Brute Force sampling", "Monte Carlo importance sampling"]

for idx, data in enumerate(dataset):
    plt.plot(
        np.log10(data[:, 0]),
        data[:, 1],
        marker=markers[idx],
        alpha=0.5,
        markersize=15,
        markeredgewidth=5,
        ls="None",
        label=labels[idx],
    )

    plt.errorbar(
        x=np.log10(data[:, 0]),
        y=data[:, 1],
        ls="None",
        yerr=data[:, 2],
        capsize=5,
        color="k",
        elinewidth=2,
    )


plt.xlabel("$Log_{10}$(x)", fontsize=28)
plt.ylabel("Value of the integral", fontsize=28)
plt.title(
    "Comparison of errors from different MC methods",
    fontsize=30,
)
plt.yticks(np.arange(0, 35, 5), np.arange(0, 35, 5), fontsize=24)
plt.xticks(np.arange(0, 9), np.arange(9), fontsize=24)
plt.legend(prop={"size": 25})

plt.show()
