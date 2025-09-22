import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# load table (no headers, whitespace separated example)
data = pd.read_csv("bladedat.txt", sep=r"\s+", header=None)

# convert to numpy array
Z = data.values   # shape (rows, cols)

# build x, y grid
ny, nx = Z.shape
x = np.arange(nx)
y = np.arange(ny)
X, Y = np.meshgrid(x, y)

# contour plot
plt.contourf(X, Y, Z, levels=20, cmap='viridis')  # filled contours
plt.colorbar(label="Value")
plt.xlabel("Column index")
plt.ylabel("Row index")
plt.title("2D Contour Plot of Table")
plt.show()
