import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file1 = "data_a.dat"
file2 = "data_n.dat"


a = np.loadtxt(file1, delimiter=" ", skiprows=1)
n = np.loadtxt(file2, delimiter=" ", skiprows=1)

plt.subplot(1, 2, 1)
plt.plot(a[:, 0], a[:, 1], ls="-.", lw=3)
plt.scatter(a[:, 0], a[:, 1], s=100, marker="x", color="r", linewidths=2)
plt.xlabel("a(n=40)", fontsize=18)
plt.ylabel("Energy", fontsize=18)
plt.grid()


plt.subplot(1, 2, 2)
plt.plot(n[:, 0], n[:, 1], ls="-.", lw=3)
plt.scatter(n[:, 0], n[:, 1], s=100, marker="x", color="r", linewidths=2)
plt.xlabel("n(a=5)", fontsize=18)
plt.ylabel("Energy", fontsize=18)

plt.suptitle(
    "Variation of energy with control parameters (a,n) for the square wave potential",
    fontsize=22,
)
plt.grid()
plt.show()
