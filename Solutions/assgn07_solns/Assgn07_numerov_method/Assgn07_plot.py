import numpy as np
import matplotlib.pyplot as plt

filename = "ss.out"
data = np.loadtxt(filename)

xax = data[:, 0]
yax = data[:, 1]

plt.plot(xax, yax)
plt.xlabel("Position $x$", fontsize=18)
plt.ylabel("$\psi(x)$", fontsize=18)
plt.title("Eigenfunction for n=2 (energy = 2.5)", fontsize=20)
plt.show()
