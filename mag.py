import magpylib as magpy
from magpylib.source.current import Circular, Line
from numpy import linspace
import numpy as np
import matplotlib.pyplot as plt
from magpylib import Collection, displaySystem
from tqdm import tqdm
from math import cos, sin, pi
from copy import deepcopy

def Coil(current, radius, windingsperlayer, layers, layerspace):
    l = [Circular(curr=current, dim=radius, pos=[0, 0, z]) for z in linspace(0, 5, windingsperlayer)]
    c = magpy.Collection(l)
    if layers > 1:
        for layer in range(1,layers):
            c.addSources([Circular(curr=current, dim=radius+layer*layerspace, pos=[0, 0, z]) for z in linspace(0, 5, N)])
    return c

def getQuadCoil(radius, angle, length, pos=(0,0,0), layers=1, turn_spacing=0.1, layer_spacing=0.1, turnperlayer=1, current=0):
    lines = []
    x0,y0,z0 = pos
    r = radius
    K = magpy.Collection()
    for layer in range(0, layers):
        points= []
        #Generate half bottom
        for dtheta in linspace(0, angle/2):
            x = r*cos(dtheta)
            y = r*sin(dtheta)
            points.append((x+x0, y+y0, z0))
        #Generate sides
        points.append((x+x0, y+y0, length+z0))

        #Generate top
        for dtheta in linspace(angle/2, -angle/2):
            x = r*cos(dtheta)
            y = r*sin(dtheta)
            points.append((x+x0, y+y0, length+z0))

        #Generate side
        points.append((x+x0, y+y0, z0))

        #Generate other half of bottom
        for dtheta in linspace(-angle/2, 0):
            x = r*cos(dtheta)
            y = r*sin(dtheta)
            points.append((x+x0, y+y0, z0))

        L = Line(current, vertices=points)
        K.addSources(L)
        angle = (2*layer_spacing + r*angle)/r
        z0= z0-layer_spacing
        length+=2*layer_spacing

    for turn in range(0,turnperlayer):
        Kt = deepcopy(K)
        Kt.move([turn * turn_spacing, 0, 0])
        K.addSources(Kt)
    return K


#Gets quadupole magnet
def getQuad(R, angle, length, layers=1, turn_spacing=.33, layer_spacing=.33, turnperlayer=1, current=0):
    K1 = getQuadCoil(R, angle, length, layers=layers, turn_spacing=turn_spacing, layer_spacing=layer_spacing, turnperlayer=turnperlayer, current= current)
    K1.rotate(180, axis=[0,0,1], anchor=[R,0,0])
    K2 = getQuadCoil(R, angle, length , layers=layers, turn_spacing=turn_spacing, layer_spacing=layer_spacing, turnperlayer=turnperlayer, current= -current)
    K2.rotate(180, axis=[0,0,1], anchor=[R,0,0])

    K3 = deepcopy(K1)
    K4 = deepcopy(K2)

    K2.rotate(90, axis=[0,0,1], anchor=[0,0,0])
    K3.rotate(180, axis=[0,0,1], anchor=[0,0,0])
    K4.rotate(270, axis=[0,0,1], anchor=[0,0,0])
    return magpy.Collection(K1, K2, K3, K4)
#Gets sextupole magnet
def getSec(R, angle, length, layers=1, turn_spacing=.33, layer_spacing=.33, turnperlayer=1, current=0):
    S1 = getQuadCoil(R, angle, length, layers=layers, turn_spacing=turn_spacing, layer_spacing=layer_spacing, turnperlayer=turnperlayer, current=current)
    N1 = getQuadCoil(R, angle, length, layers=layers, turn_spacing=turn_spacing, layer_spacing=layer_spacing, turnperlayer=turnperlayer, current=-current)

    S2 = deepcopy(S1)
    S3 = deepcopy(S1)
    N2 = deepcopy(N1)
    N3 = deepcopy(N1)

    S1.rotate(30, axis=[0,0,1])
    S2.rotate(150, axis=[0,0,1])
    S3.rotate(270, axis=[0,0,1])

    N1.rotate(90, axis=[0,0,1])
    N2.rotate(210, axis=[0,0,1])
    N3.rotate(330, axis=[0,0,1] )
    return magpy.Collection(N1, N2, N3, S1, S2, S3)

#mm, mT
R = 10
N = 50
layer_spacing = .33 #mm
target_angle = pi/20
angle = (R*target_angle - (2*N*layer_spacing))/R
length = 10
layers = 50
current = .00001

# N1 = Circular(curr=current*layers*N, dim=R)
# N1.rotate(90, axis=[1,0,0])
# N1.move([0,-R,0])
# N2 = deepcopy(N1)
#
# S1 = Circular(curr=-current*layers*N, dim=R)
# S1.rotate(90, axis=[1,0,0])
# S1.move([0,-R,0])
# S2 = deepcopy(S1)
#
# N2.rotate(180,axis=[0,0,1], anchor=[0,0,0])
# S1.rotate(90, axis=[0,0,1], anchor=[0,0,0])
# S2.rotate(270,axis=[0,0,1], anchor=[0,0,0])

#c = magpy.Collection(N1, N2, S1, S2)
c = getQuad(R, pi/6, length, layers=1, turn_spacing=.33, layer_spacing=layer_spacing, turnperlayer=1, current=.01)
#c = getQuad(R, angle, length, layers=1, turn_spacing=.33, layer_spacing=layer_spacing, turnperlayer=1, current=1*N*layers)

#c = getSec(R, angle, length, layers=1, turn_spacing=.33, layer_spacing=layer_spacing, turnperlayer=1, current=10)
#g = getSec(R, angle, length, layers=1, turn_spacing=.33, layer_spacing=layer_spacing, turnperlayer=1, current=1*N*layers)
#g.move([0,0,-(200+(2*layer_spacing*100))])
#c = magpy.Collection(c,g)

print("Source Generation Complete")

# import pickle
# f = open('quad', 'rb')
# c = pickle.load(f)
# f.close()
#


# calculate B-field on a grid
points = 100
xs = np.linspace(-R,R,points)
zs = np.linspace(-R,R,points)
POS = np.array([(x,z,0) for z in zs for x in xs])
print("Computing Magnetic Field...")
print(POS.size)
Bs = c.getB(POS).reshape(points,points,3)*1000     #<--VECTORIZED
print("Magnetic Field computed")

# create figure
fig = plt.figure(figsize=(12,6))
ax1 = fig.add_subplot(121, projection='3d')  # 3D-axis
ax2 = fig.add_subplot(122)                   # 2D-axis

# display system geometry on ax1
displaySystem(c, subplotAx=ax1, suppress=True)
#
# # display field in xz-plane using matplotlib
X,Z = np.meshgrid(xs,zs)
U,V = Bs[:,:,0], Bs[:,:,1]
ax2.streamplot(X, Z, U, V, color=np.log(U**2+V**2), density=2)

plt.show()

q = -1.602e-19
m = 9.109e-31
Vmag = 100000 #Voltage of accelerator plate

import math
# E = gamma*m*c^2 = mc^2/sqrt(1-B^2)
KE = 7 # Fermi energy, eV
cs = 3e8 #m/s
vmax = math.sqrt(abs((KE*q)**2 - (m*(cs**2))**2)) * cs * 1e3 #mm/s
print(vmax)

import random
particles = []
for i in tqdm(range(0,1)):
    v = np.asarray( [ random.random(), random.random(), -random.random() ] )
    v = (vmax/math.sqrt( v[0]**2 + v[1]**2 + v[2]**2 ))*v
    v[2] = -1.6e6 *1e6
    V[0] = 0
    v[1] = 0
    x = [random.randint(-1000,1000)/1000, random.randint(-1000,1000)/1000, length] #mm
    particles.append([x, v])

#
time_max = 1e-8
time_step = 1000000
for particle in tqdm(particles):
    x = particle[0]
    v = particle[1]
    V = []
    X = []
    A = []
    T = []
    V.append(v)
    X.append(x)
    tv = linspace(0, time_max, time_step)  # s
    h = time_max / time_step
    z = x[2]
    i = 0
    while z > -10 and z <= particle[0][2] + 1e-6 and abs(x[0]) < 3*R and abs(x[1]) < 3*R:
        z = x[2]
        try:
            t = tv[i]

            B = c.getB(x) * 1000  # mT

            x = x + v * h

            t = t + h
            #Runge-kutta
            k1 = (q / m) * np.cross(v, B)
            k2 = (q / m) * np.cross(v+ (h*k1)/2, B)
            k3 = (q / m) * np.cross(v+ (h*k2)/2, B)
            k4 = (q / m) * np.cross(v+ h*k3/2, B)

            v = v + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4)

            V.append(v)
            X.append(x)
            T.append(t)
            if i % 1000 == 0:
                print(f'\nPosition {x}, Velocity {v}')

            i += 1
        except(IndexError):
            print('Maxed Time Vector...')
            break

    X = np.asarray(X)

    ax1.plot3D(X[:, 0], X[:, 1], X[:, 2])
ax1.set_xlim([-30, 30])
ax1.set_ylim([-30, 30])