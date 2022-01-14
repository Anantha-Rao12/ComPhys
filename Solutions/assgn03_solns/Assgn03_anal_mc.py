import numpy as np
import matplotlib.pyplot as plt

print(
    "Hello!\n This program reads and analyses the data is file named as 'ising_3d_T_L20_random.dat'"
)

data = np.loadtxt("ising_3d_T_L20_random.dat")

print(
    "Magnetization per spin:", np.mean(data, axis=0)[1], "+-", np.std(data, axis=0)[1]
)
print("Energy per spin:", np.mean(data, axis=0)[2], "+-", np.std(data, axis=0)[-1])

# plt.plot(data.T[0], data.T[1], "r--", alpha=0.4)
# plt.show()
