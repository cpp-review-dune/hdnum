import matplotlib.pyplot as plt 
import numpy as np

iteration_points = []

Filename1 = "projectedNewtonSolver.dat" # result of projected newton method
Filename2 = "loss.dat" # loss(x) = f(x)
Filename3 = "constraints.dat" # constraints of type lower bound <= Ax <= upper bound

f = open(Filename1,"r")

# get the iteration points
for i,a in enumerate(f):
    x = a.split("   ")
    iteration_points.append(np.array([np.double(x[0]),np.double(x[1])]))
        
iteration_points = np.array(iteration_points, dtype= np.double)

f = open(Filename2, "r")

losses = []
minValue = 0
maxValue = 0

# get the losses
for i,a in enumerate(f):
    if i == 0: # the first line contains the domain
        x = a.split("   ")
        minValue = np.double(x[0])
        maxValue = np.double(x[1])
    else:
        losses.append(np.double(a))

domainData = np.arange(minValue, maxValue, 0.2)
dimension = int(np.sqrt(len(losses)))

if len(domainData) != dimension: # fix bug issue
    domainData = np.arange(minValue, maxValue + 0.0001, 0.2)

X,Y = np.meshgrid(domainData,domainData)
Z = np.reshape(losses, (dimension, dimension)).T


f = open(Filename3, "r")
lowerbounds = []
upperbounds = []
constraints = []

# get constraints
for i,a in enumerate(f):
    x = a.split("   ")
    lowerbounds.append(np.double(x[0]))
    constraints.append(np.array([np.double(x[1]),np.double(x[2])]))
    upperbounds.append(np.double(x[3]))

constraints = np.array(constraints, dtype=np.double)

# compute the feasible space
feasibleSpace = Z == Z 
for i in range(constraints.shape[0]):
    constrainFulfilled = np.logical_and(constraints[i][0] * X + constraints[i][1] * Y >= lowerbounds[i], constraints[i][0] * X + constraints[i][1] * Y <= upperbounds[i])
    feasibleSpace = np.logical_and(feasibleSpace, constrainFulfilled)

fig, ax = plt.subplots(figsize= (9,9))
ax.contourf(X,Y,feasibleSpace, alpha = 0.2)
ax.contour(X,Y,Z, 200, alpha=0.2)
ax.scatter(iteration_points[:,0], iteration_points[:,1])
ax.scatter(iteration_points[0][0], iteration_points[0][1], color="r")
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.show()
