import numpy as np
import matplotlib.pyplot as plt

Filename = "loss.dat"
Filename1= "newton_solver.dat"
Filename2= "dog_leg_solver.dat"

#filename
f = open(Filename1, "r")
iteration_points = []
directions = []
step_sizes = []
norms= []
reductions = []
loss = []


for i,a in enumerate(f):
    if i >= 2:
        x = a.split("   ")
        iteration_points.append(np.array([np.double(x[0]),np.double(x[1])]))
        directions.append(np.array([np.double(x[2]),np.double(x[3])]))
        step_sizes.append(np.double(x[4]))
        norms.append(np.double(x[5]))
        reductions.append(np.double(x[6]))
        loss.append(np.double(x[7]))
        
iteration_points = np.array(iteration_points, dtype= np.double)
directions = np.array(directions, dtype= np.double)

f = open(Filename2, "r")

iteration_points_dog_leg = []
newton_directions_dog_leg = []
steepest_descent_directions_dog_leg = []
dog_leg_directions_dog_leg = []
trust_radii = []
norms_dog_leg= []
reductions_dog_leg = []
loss_dog_leg = []

for i,a in enumerate(f):
    if i >= 2:
        x = a.split("   ")
        iteration_points_dog_leg.append(np.array([np.double(x[0]),np.double(x[1])]))
        newton_directions_dog_leg.append(np.array([np.double(x[2]),np.double(x[3])]))
        steepest_descent_directions_dog_leg.append(np.array([np.double(x[4]),np.double(x[5])]))
        dog_leg_directions_dog_leg.append(np.array([np.double(x[6]),np.double(x[7])]))
        trust_radii.append(np.double(x[8]))
        norms_dog_leg.append(np.double(x[9]))
        reductions_dog_leg.append(np.double(x[10]))
        loss_dog_leg.append(np.double(x[11]))

iteration_points_dog_leg = np.array(iteration_points_dog_leg, dtype = np.double)
newton_directions_dog_leg = np.array(newton_directions_dog_leg, dtype = np.double)
steepest_descent_directions_dog_leg = np.array(steepest_descent_directions_dog_leg, dtype = np.double)
dog_leg_directions_dog_leg = np.array(dog_leg_directions_dog_leg, dtype = np.double)

f = open(Filename, "r")
min_value = 0
max_value = 0
res = 0
list = []

for i,a in enumerate(f):
    if i == 0:
        x = a.split(" ")
        min_value = np.double(x[0])
        max_value = np.double(x[1])
    else:
        list.append(np.double(a))

l = np.arange(min_value, max_value, 0.2)
dimension = int(np.sqrt(len(list)))

if len(l) != dimension:
    l = np.arange(min_value, max_value + 0.0001, 0.2)

X,Y = np.meshgrid(l,l)
Z = np.reshape(list, (dimension, dimension)).T


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')
surf = ax.plot_wireframe(X, Y, Z)
ax.scatter(iteration_points[:,0], iteration_points[:,1], loss,  color = "k", s=10)
ax.scatter(iteration_points_dog_leg[:,0], iteration_points_dog_leg[:,1], loss_dog_leg,  color = "r", s=10)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.legend(["Surface of loss function" ,"Newton", "Newton dog leg cauchy"])
plt.show()

overlapping = 0.1
plt.contour(X,Y,Z,200, alpha=overlapping)
plt.scatter(iteration_points[:,0], iteration_points[:,1], color = "b")
plt.scatter(iteration_points_dog_leg[:,0], iteration_points_dog_leg[:,1], color = "orange")
plt.scatter(iteration_points[0,0], iteration_points[0,1], color="red")
plt.xlabel("x")
plt.ylabel("y")
plt.legend(['Newton', 'Newton dog leg cauchy',"Initial solution"])
plt.show()

plt.plot(norms)
plt.plot(norms_dog_leg)
plt.xlabel("number of iterations")
plt.ylabel("norm of residial")
plt.legend(['Newton', 'Newton dog leg cauchy'])
plt.show()

plt.plot(reductions)
plt.plot(reductions_dog_leg)
plt.xlabel("number of iterations")
plt.ylabel("reduction")
plt.legend(['Newton', 'Newton dog leg cauchy'])
plt.show()

fig, ax = plt.subplots()
overlapping = 0.1
ax.contour(X,Y,Z,200, alpha=overlapping)

for i in range(dog_leg_directions_dog_leg.shape[0]-1):
    x = iteration_points_dog_leg[i]
    newton = newton_directions_dog_leg[i]
    steepest = steepest_descent_directions_dog_leg[i]
    dog_leg = dog_leg_directions_dog_leg[i]

    trust_radius = trust_radii[i]

    x_newton = x + newton
    x_steepest = x + steepest
    x_dog_leg = x + dog_leg

    circle1 = plt.Circle((x[0], x[1]), trust_radius, color='k', fill=False)
    ax.add_patch(circle1)
    ax.set_aspect('equal')

    ax.scatter(x[0], x[1], color = "k")

    """
    #to vizualize the dog leg direction
    if dog_leg[0] != 0 and dog_leg[1] != 0:
        ax.plot(np.array([x[0], x_newton[0]]), np.array([x[1], x_newton[1]]), color ="b")
        ax.plot(np.array([x[0], x_steepest[0]]), np.array([x[1], x_steepest[1]]) , color= "r")
        ax.plot(np.array([x[0], x_dog_leg[0]]), np.array([x[1], x_dog_leg[1]]), color = "g" )
        ax.plot(np.array([x_steepest[0], x_steepest[0] + (x_newton-x_steepest)[0]]), np.array([x_steepest[1], x_steepest[1] + (x_newton-x_steepest)[1]]))
        break
    """

    ax.plot(np.array([x[0], x_newton[0]]), np.array([x[1], x_newton[1]]), color ="b")
    ax.plot(np.array([x[0], x_steepest[0]]), np.array([x[1], x_steepest[1]]) , color= "r")
    ax.plot(np.array([x[0], x_dog_leg[0]]), np.array([x[1], x_dog_leg[1]]), color = "g" )

plt.xlabel("x")
plt.ylabel("y")
plt.legend(["trust radius", "current solution", "newton direction", "steepest descent direction", "dog leg direction"])
plt.show()