import numpy as np
import matplotlib.pyplot as plt

files = "q3_sinx_x0_v-1e-1.dat q3_sinx_x0_v3e-2.dat q3_sinx_x0_v3e-1.dat q3_x_x0_v-1e-1.dat q3_x_x0_v3e-2.dat q3_x_x0_v3e-1.dat".split()

titles = []
for f in files:
    vals = f.split("_")
    title = (
        r"$d^2x/dt^2 =$ " + f"{vals[1]}; $x_0$={vals[2][1:]}; $v_0$={vals[-1][1:-4]} "
    )
    titles.append(title)

labels = [r"$d^2x/dt^2 = -sinx$"] * 3 + [r"$d^2x/dt^2 = -x$"] * 3

for i in range(6):
    data = np.loadtxt(files[i])
    plt.subplot(2, 3, i % 3 + 1)
    plt.plot(data[:, 0], data[:, 1], label=labels[i])
    plt.xlabel("Time (t)", fontsize=18)
    plt.ylabel("Position (x)", fontsize=18)
    plt.legend()

    print(i % 3 + 4)
    plt.subplot(2, 3, i % 3 + 4)
    plt.plot(data[:, 0], data[:, 2], label=labels[i])
    plt.xlabel("Time (t)", fontsize=18)
    plt.ylabel("Velocity (v)", fontsize=18)
    plt.legend()

plt.show()


"""for idx, dataset in enumerate(files):
    plt.subplot(2, 3, idx + 1)
    print(dataset)
    data = np.loadtxt(dataset)
    plt.plot(data[:, 1], data[:, 2], label="Trajectory", color="royalblue")
    plt.plot(
        data[0, 1],
        data[0, 2],
        marker="x",
        markersize=10,
        color="limegreen",
        markeredgewidth=3,
        label="Starting point",
    )
    plt.xlabel("Position (x)", fontsize=18)
    plt.ylabel("Velocity (v)", fontsize=18)

    # plt.plot(data[:, 0], data[:, 2], label="v")
    plt.legend()
    plt.title(titles[idx])

plt.show()
"""
