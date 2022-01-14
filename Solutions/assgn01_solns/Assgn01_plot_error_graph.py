import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

data = np.loadtxt("trap_pi.dat")

want = data[:9, (0, -1)]

p = np.polyfit(want[:, 0], want[:, -1], 1)
print(p)


r2 = r2_score(want[:, -1], p[0] + want[:, 0] * p[1])
plt.plot(
    want[:, 0],
    p[0] + want[:, 0] * p[1],
    "k--",
    label="Fitted Line \n" + f"$y={p[1]:.3}x{p[0]:.3}$\n" + f"$R^2$:{r2:.3}",
)

plt.plot(data[:, 0], data[:, -1], "b--")
plt.plot(data[:, 0], data[:, -1], "rp", markersize=7)
plt.xlabel("Log(No of bins)", fontsize=18)
plt.ylabel("Log(Abs(Error))", fontsize=18)
plt.title(
    "Composite Trapezoidal integration for " + r"$\int_0^1 \frac{4}{(1+x^2)}$",
    fontsize=20,
)

plt.xticks(np.arange(2.5, 21, 2.5), fontsize=14)
plt.yticks(np.arange(-5, -31, -5), fontsize=14)

plt.legend(prop={"size": 15})
plt.show()
