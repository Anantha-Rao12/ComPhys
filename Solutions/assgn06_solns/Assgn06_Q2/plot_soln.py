import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("b_value_dx_01_limit_0001.dat")

xax = [f"{i:.2}" for i in np.arange(0, 1.1, 0.2)]
yax = [int(i) for i in np.arange(0, 41.1, 5)]

plt.scatter(data[:, 0], data[:, 1], marker="x")
plt.xlabel("X", fontsize=20)
plt.ylabel("Y(x)", fontsize=20)
plt.title("Q2. Solution of Diff Eq: with Boundary conditions with dx=0.01", fontsize=20)

plt.xticks(np.arange(0, 1.1, 0.2), xax, fontsize=18)
plt.yticks(np.arange(0, 41.1, 5), yax, fontsize=18)
plt.grid()
plt.show()
