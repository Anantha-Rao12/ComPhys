! Code to Analyse the results of the MC simulation to solve the 3D Isiing model
! Author: Anantha RAo (Reg no:20181044)


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_df(path: str) -> pd.DataFrame:
"""Reads the data file from the simulation using the file path and return a pandas dataframe object"""
    data = np.loadtxt(path)
    return pd.DataFrame(data, columns=["T", "<M>", "<E>", "cv", "chi", "B_c"])


df7 = get_df("ising_L7_3d_data_new.dat")
df8 = get_df("ising_L8_3d_data_new.dat")
df9 = get_df("ising_L9_3d_data_new.dat")
df = [df7, df8, df9]
nspins = [7, 8, 9]
y_label = [
    "Magnetization per spin " + r"$\langle \frac{M}{N} \rangle$",
    "Energy per spin " + r"$\langle \frac{E}{N} \rangle$",
    "$C_v$(T)",
    r"$\chi (T)$",
]
markers = ["o", "s", "P"]

# print(df7.tail())
# want = np.where(df7.T.values == 3.8)[0]
# print(want)
# print(df7.iloc[want, :])

"""  Comment to plot the simulation data
fig, axes = plt.subplots(2, 2, dpi=100)

for idx, ax in enumerate(axes.reshape(-1)):
    for n_spins in range(3):
        ax.plot(
            df[n_spins].iloc[:, 0],
            df[n_spins].iloc[:, idx + 1],
            label=f"L={nspins[n_spins]}",
            lw=3,
            marker=markers[n_spins],
            markersize=5,
            alpha=0.5,
        )
    ax.legend(prop={"size": 20})
    ax.grid()
    ax.set_ylabel(y_label[idx], fontsize=16)
    ax.set_xlabel("Temperature (T) (kB units)", fontsize=16)
    ax.set_xticks(np.arange(3.8, 4.75, 0.15))
    ax.set_xticklabels([f"{i:.2}" for i in np.arange(3.8, 4.75, 0.15)], fontsize=14)

plt.suptitle("3D Ising Model Simulation Results", fontsize=24)
plt.show()


"""

""" Comment to plot the Binders cumulant data
# Plot Binders Cumulant

for i in range(3):
    plt.plot(
        df[i].iloc[:, 0],
        df[i].iloc[:, -1],
        label=f"L={nspins[i]}",
        lw=5,
        alpha=0.6,
        marker="o",
    )
plt.axvline(x=4.52, label="x=4.52", lw=3, alpha=0.5)
plt.axvline(x=4.50, label="x=4.50", lw=5, alpha=0.5, color="k")
plt.axvline(x=4.48, label="x=4.48", lw=3, alpha=0.5)

plt.xlabel("Temperature (T) in kB units", fontsize=18)
plt.ylabel(
    "$U_L$" + r"=$1 - \frac{\langle M^4 \rangle _L}{3 \langle M^2 \rangle ^2_L}$",
    fontsize=20,
)

plt.title("Finding Critical Temperature from Binder's Cumulant", fontsize=24)
plt.legend(prop={"size": 20})
plt.xticks(np.arange(3.8, 4.76, 0.075), fontsize=18)
plt.yticks(fontsize=18)
plt.grid()
plt.show()


"""
