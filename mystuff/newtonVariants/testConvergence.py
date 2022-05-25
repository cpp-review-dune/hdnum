Filename = "convergence.dat"  # contains results of type: has converged or not
Filename2 = "convergence2.dat" # contains results of type: to which optimum did the method converge

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

f = open(Filename, "r")

iterationsNeeded = [] # number of iterations needed to converge. -1 means no convergence.

fig, ax = plt.subplots(1,2, figsize=(15,8))

minRange = 0
maxRange = 0

mode = 0 # 0 stands for test for different initial solutions, 1 for different initial and maximum trust radii

for i,a in enumerate(f): # iterate over lines in file. First line contains the mode. Second line contains the range.
    if i == 0:
        mode = int(a)
    elif i == 1:
        x = a.split(" ")
        minRange = np.double(x[0])
        maxRange = np.double(x[1])
    else:
        iterationsNeeded.append(np.double(a))

if mode == 1:
    ax[0].set_xlabel("x_coordinate")
    ax[0].set_ylabel("y_coordinate")
    ax[1].set_xlabel("x_coordinate")
    ax[1].set_ylabel("y_coordinate")
else:
    ax[0].set_xlabel("Maximum trust radius")
    ax[0].set_ylabel("Initial trust radius")
    ax[1].set_xlabel("Maximum trust radius")
    ax[1].set_ylabel("Initial trust radius")

ax[0].set_title("Number of iterations to converge")

dimension = int(np.sqrt(len(iterationsNeeded))) # reshape list to quadratic matrix for visualization
iterationsNeeded = np.reshape(iterationsNeeded, (dimension, dimension)).T

im = ax[0].imshow(iterationsNeeded,interpolation='none', extent=[minRange,maxRange,maxRange, minRange])
divider = make_axes_locatable(ax[0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

f = open(Filename2, "r")

solutionIndex = [] # index of the set of all found solution in the domain. -1 means no convergence

for i,a in enumerate(f): # iterate over lines in file
    solutionIndex.append(np.double(a))

solutionIndex = np.array(solutionIndex)

dimension = int(np.sqrt(len(solutionIndex)))

solutionIndex = np.reshape(solutionIndex, (dimension, dimension)).T

ax[1].imshow(solutionIndex,interpolation='none', extent=[minRange,maxRange,maxRange, minRange])
ax[1].set_title("Solution index")
plt.show()

