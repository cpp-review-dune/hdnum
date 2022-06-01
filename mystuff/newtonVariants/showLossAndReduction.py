import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

FilenameNewton = "newton_solver.dat"
FilenameDogLeg = "dog_leg_solver.dat"

f = open(FilenameNewton, "r")

lossNewton = []
reductionNewton = []

for i,a in enumerate(f):
    if i >= 2:
        x = a.split("   ")
        lossNewton.append(np.double(x[3]))
        reductionNewton.append(np.double(x[4]))

f = open(FilenameDogLeg, "r")

lossDogLeg = []
reductionDogLeg = []

for i,a in enumerate(f):
    if i >= 2:
        x = a.split("   ")
        lossDogLeg.append(np.double(x[5]))
        reductionDogLeg.append(np.double(x[6]))

fig, ax = plt.subplots(1,2, figsize=(15,8))

ax[0].plot(lossNewton)
ax[0].plot(lossDogLeg)
ax[0].set_xlabel("number of iterations")
ax[0].set_ylabel("norm")
ax[0].legend(["Newton method", "Newton dogleg cauchy method"])

ax[1].plot(reductionNewton)
ax[1].plot(reductionDogLeg)
ax[1].set_xlabel("number of iterations")
ax[1].set_ylabel("reduction")
ax[1].legend(["Newton method", "Newton dogleg cauchy method"])

plt.show()




