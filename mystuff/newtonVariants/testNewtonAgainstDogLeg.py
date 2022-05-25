import numpy as np
import matplotlib.pyplot as plt

Filename = "loss.dat" # loss(x) = 0.5 * F(x)^T * F(x) where F(x) = 0 defines the nonlinear problem
Filename1= "newton_solver.dat" # result of newton solver
Filename2= "dog_leg_solver.dat" # result of newton dogleg cauchy solver

f = open(Filename1, "r")

iterationPointsNewton = []
directionsNewton = []
stepSizesNewton = []
normsNewton= []
reductionsNewton = []
lossesNewton = []


for i,a in enumerate(f): # iterate over lines in file. The first two lines contains information about the file.
    if i >= 2:
        x = a.split("   ")
        iterationPointsNewton.append(np.array([np.double(x[0]),np.double(x[1])]))
        directionsNewton.append(np.array([np.double(x[2]),np.double(x[3])]))
        stepSizesNewton.append(np.double(x[4]))
        normsNewton.append(np.double(x[5]))
        reductionsNewton.append(np.double(x[6]))
        lossesNewton.append(np.double(x[7]))
        
iterationPointsNewton = np.array(iterationPointsNewton, dtype= np.double)
directionsNewton = np.array(directionsNewton, dtype= np.double)

f = open(Filename2, "r")

iterationPointsDogleg = []
newtonDirectionsDogleg = []
steepestDescentDirections = []
doglegDirections = []
trustRadii = []
normsDogleg= []
reductionsDogleg = []
lossesDogleg = []

for i,a in enumerate(f): # iterate over lines in file. The first two lines contains information about the file.
    if i >= 2:
        x = a.split("   ")
        iterationPointsDogleg.append(np.array([np.double(x[0]),np.double(x[1])]))
        newtonDirectionsDogleg.append(np.array([np.double(x[2]),np.double(x[3])]))
        steepestDescentDirections.append(np.array([np.double(x[4]),np.double(x[5])]))
        doglegDirections.append(np.array([np.double(x[6]),np.double(x[7])]))
        trustRadii.append(np.double(x[8]))
        normsDogleg.append(np.double(x[9]))
        reductionsDogleg.append(np.double(x[10]))
        lossesDogleg.append(np.double(x[11]))

iterationPointsDogleg = np.array(iterationPointsDogleg, dtype = np.double)
newtonDirectionsDogleg = np.array(newtonDirectionsDogleg, dtype = np.double)
steepestDescentDirections = np.array(steepestDescentDirections, dtype = np.double)
doglegDirections = np.array(doglegDirections, dtype = np.double)

f = open(Filename, "r")
minRange = 0
maxRange = 0
res = 0
losses = []

for i,a in enumerate(f): # iterate over lines in file. The first line contains the range.
    if i == 0:
        x = a.split(" ")
        minRange = np.double(x[0])
        maxRange = np.double(x[1])
    else:
        losses.append(np.double(a))

domainData = np.arange(minRange, maxRange, 0.2)
dimension = int(np.sqrt(len(losses)))

if len(domainData) != dimension: # fix bug issue
    domainData = np.arange(minRange, maxRange + 0.0001, 0.2)

X,Y = np.meshgrid(domainData,domainData)
Z = np.reshape(losses, (dimension, dimension)).T 


# 3D visualization
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')
surf = ax.plot_wireframe(X, Y, Z)
ax.scatter(iterationPointsNewton[:,0], iterationPointsNewton[:,1], lossesNewton,  color = "k", s=10)
ax.scatter(iterationPointsDogleg[:,0], iterationPointsDogleg[:,1], lossesDogleg,  color = "r", s=10)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.legend(["Surface of loss function" ,"Newton", "Newton dog leg cauchy"])
plt.show()

# 2D visualization using contourplots
overlapping = 0.1
plt.contour(X,Y,Z,200, alpha=overlapping)
plt.scatter(iterationPointsNewton[:,0], iterationPointsNewton[:,1], color = "b")
plt.scatter(iterationPointsDogleg[:,0], iterationPointsDogleg[:,1], color = "orange")
plt.scatter(iterationPointsNewton[0,0], iterationPointsNewton[0,1], color="red") # initial solution
plt.xlabel("x")
plt.ylabel("y")
plt.legend(['Newton', 'Newton dog leg cauchy',"Initial solution"])
plt.show()

# trajectory of the norm
plt.plot(normsNewton)
plt.plot(normsDogleg)
plt.xlabel("number of iterations")
plt.ylabel("norm of residial")
plt.legend(['Newton', 'Newton dog leg cauchy'])
plt.show()

# trajectory of the reduction
plt.plot(reductionsNewton)
plt.plot(reductionsDogleg)
plt.xlabel("number of iterations")
plt.ylabel("reduction")
plt.legend(['Newton', 'Newton dog leg cauchy'])
plt.show()

fig, ax = plt.subplots()
overlapping = 0.1
ax.contour(X,Y,Z,200, alpha=overlapping)

# Vizualize the trust regions
for i in range(doglegDirections.shape[0]-1):
    x = iterationPointsDogleg[i]
    newton = newtonDirectionsDogleg[i]
    steepest = steepestDescentDirections[i]
    dogleg = doglegDirections[i]

    trustRadius = trustRadii[i]

    xNewton = x + newton
    xSteepest = x + steepest
    xDogleg = x + dogleg

    circle1 = plt.Circle((x[0], x[1]), trustRadius, color='k', fill=False)
    ax.add_patch(circle1)
    ax.set_aspect('equal')

    ax.scatter(x[0], x[1], color = "k")

    ax.plot(np.array([x[0], xNewton[0]]), np.array([x[1], xNewton[1]]), color ="b")
    ax.plot(np.array([x[0], xSteepest[0]]), np.array([x[1], xSteepest[1]]) , color= "r")
    ax.plot(np.array([x[0], xDogleg[0]]), np.array([x[1], xDogleg[1]]), color = "g" )

plt.xlabel("x")
plt.ylabel("y")
plt.legend(["trust radius", "current solution", "newton direction", "steepest descent direction", "dog leg direction"])
plt.show()
