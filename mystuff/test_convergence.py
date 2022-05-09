Filename = "convergence.dat"

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
f = open(Filename, "r")

list = []

fig, ax = plt.subplots()

min_range = 0
max_range = 0

mode = 0
for i,a in enumerate(f):
    if i == 0:
        mode = int(a)
    elif i == 1:
        x = a.split(" ")
        min_range = np.double(x[0])
        max_range = np.double(x[1])
    else:
        list.append(np.double(a))

if mode == 1:
    ax.set_xlabel("x_coordinate")
    ax.set_ylabel("y_coordinate")
else:
    ax.set_xlabel("Maximum trust radius")
    ax.set_ylabel("Initial trust radius")

ax.set_title("Number of iterations to converge")


dimension = int(np.sqrt(len(list)))
list = np.reshape(list, (dimension, dimension)).T

im = ax.imshow(list,interpolation='none', extent=[min_range,max_range,max_range, min_range])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

plt.show()
