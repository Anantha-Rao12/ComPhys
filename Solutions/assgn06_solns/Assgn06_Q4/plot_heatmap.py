import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data = np.loadtxt("data.dat")
temp = np.zeros((34, 34))
df = pd.DataFrame(data, columns=["x", "y", "t"])

for index, row in df.iterrows():
    print(row.x, row.y, row.t)
    temp[int(row.x) - 1, int(row.y) - 1] = row.t
    # temp[row.x - 1, row.y - 1] = row.t

sns.heatmap(temp.T, annot=False)
plt.xlabel("X", fontsize=18)
plt.ylabel("Y", fontsize=18)
plt.title("Temperature profile for the 2D lattice (Q4)", fontsize=20)
plt.xticks(np.arange(0.5, 34.5, 5), np.arange(1, 35, 5))
plt.yticks(np.arange(0.5, 34.5, 5), np.arange(1, 35, 5))
plt.show()
