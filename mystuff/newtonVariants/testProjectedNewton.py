import matplotlib.pyplot as plt 
import numpy as np

iteration_points = []
file = "a.dat"
file2 = "loss.dat"
file3 = "constraints.dat"

f = open(file, "r")

for i,a in enumerate(f):
    x = a.split("   ")
    iteration_points.append(np.array([np.double(x[0]),np.double(x[1])]))
        
iteration_points = np.array(iteration_points, dtype= np.double)

f = open(file2, "r")

list = []

for i,a in enumerate(f):
    list.append(np.double(a))

l = np.arange(-5, 10, 0.2)
dimension = int(np.sqrt(len(list)))

if len(l) != dimension:
    l = np.arange(-5, 10 + 0.0001, 0.2)

X,Y = np.meshgrid(l,l)
Z = np.reshape(list, (dimension, dimension)).T


f = open(file3, "r")
lowerbounds = []
upperbounds = []
constraints = []
for i,a in enumerate(f):
    x = a.split("   ")
    lowerbounds.append(np.double(x[0]))
    constraints.append(np.array([np.double(x[1]),np.double(x[2])]))
    upperbounds.append(np.double(x[3]))

constraints = np.array(constraints, dtype=np.double)
#print(lowerbounds)
#print(upperbounds)
#print(constraints)

P = Z == Z
#print(P)
for i in range(constraints.shape[0]):
    O = np.logical_and(constraints[i][0] * X + constraints[i][1] * Y >= lowerbounds[i], constraints[i][0] * X + constraints[i][1] * Y <= upperbounds[i])
    P = np.logical_and(P,O)

#print(P)

plt.contourf(X,Y,P, alpha = 0.2)
plt.contour(X,Y,Z,200, alpha=0.2)
plt.scatter(iteration_points[:,0], iteration_points[:,1])
plt.show()
