import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("D_erfvalues.dat")

plt.plot(data[:, 0], data[:, 1], "k--", lw=5)
plt.grid()
plt.xlabel("x", fontsize=16)
plt.ylabel("Erf(x)", fontsize=16)
plt.title(
    "Plotting error function from trapezoidal integration of "
    + r"$\frac{2}{\sqrt{\pi}}\int_0^x e^{-t^2}dt$",
    fontsize=20,
)
plt.show()
