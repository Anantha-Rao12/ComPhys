import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Problem1 n=100 data

print("Reading data for 100 RNs case")
p11_randompts = np.loadtxt("p1_randomnumbers_100.dat")
p11_ck = np.loadtxt("p1_ckvalues_100.dat")
p11_scatter = np.loadtxt("p1_scattervalues_100.dat")


# Problem1 n=10000 data

print("Reading data for 10000 RNs case")
p12_randompts = np.loadtxt("p1_randomnumbers_10_000.dat")
p12_ck = np.loadtxt("p1_ckvalues_10_000.dat")
p12_scatter = np.loadtxt("p1_scattervalues_10_000.dat")

# Problem1 n=10_00_000 data

print("Reading data for 1000000 RNs case")
p13_randompts = np.loadtxt("p1_randomnumbers_10_00_000.dat")
p13_ck = np.loadtxt("p1_ckvalues_10_00_000.dat")
p13_scatter = np.loadtxt("p1_scattervalues_10_00_000.dat")


# Plots

hist = [
    p11_randompts,
    p12_randompts,
    p13_randompts,
]

mean = list(map(np.mean, hist))
std = list(map(np.std, hist))

scatter = [
    p11_scatter,
    p12_scatter,
    p13_scatter,
]

ck = [p11_ck, p12_ck, p13_ck]
npts = [100, 10000, 1000000]


fig, axes = plt.subplots(nrows=3, ncols=3)
for idx, ax in enumerate(axes):
    ax[0].hist(hist[idx][:], rwidth=0.92, bins=20)
    ax[0].text(
        0.5,
        0.65,
        f"Mean = {mean[idx]:.4} \n Std.dev = {std[idx]:.4}",
        transform=ax[0].transAxes,
    )
    ax[0].set_title(f"Histogram of {npts[idx]} random numbers", fontsize=14)

    ax[1].plot(scatter[idx][:, 0], scatter[idx][:, 1], "ro", alpha=0.5)
    ax[1].set_title("Scatter plot", fontsize=14)
    ax[1].set_xlabel(r"$X_{i}$", fontsize=14)
    ax[1].set_ylabel(r"$X_{i+1}$", fontsize=14)

    ax[2].plot(ck[idx][:, 0], ck[idx][:, 1])
    ax[2].set_title("Autocorrelation ($C_k$) Plot", fontsize=14)
    ax[2].set_xlabel("k", fontsize=14)
    ax[2].set_ylabel("Correlation value", fontsize=14)


plt.tight_layout()
plt.show()
