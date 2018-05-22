import astropy
import numpy as np
import matplotlib.pyplot as plt

global G
global m_sun
global m_jup
G = 6.674e-11 #m^3 kg ^-1 s ^-2
m_sun = 1.989e30 #kg
m_jup = 1.898e27 #kg
dmax = 8.157e11 #max distance from sun
dmin = 7.409e11 #min distance from jupiter to sun
v_jup = 12.44e3 #ms^-1


semimajor = 0.5*(dmax+dmin) #length of semimajor axis in m
T = 2*np.pi*np.sqrt((semimajor**3)/(G*m_sun)) #calculated period of orbit 
tmax = 6*T #seconds

nsteps = 50000
nbodies = 3
dt = tmax/nsteps

#generate empty arrays for x, y, z
x = np.empty([nsteps, 3*nbodies])
v = np.empty([nsteps, 3*nbodies])
a = np.empty([nsteps, 3*nbodies])
m = np.empty([nbodies])

#initial values for x, v
x[0] = [0.0, 0.0, 0.0, dmax, 0.0, 0.0, 0.5*dmax, np.sqrt(3)*0.5*dmax, 0.0]
v[0] = [0.0, 0.0, 0.0, 0.0, v_jup, 0.0, 0.0, -np.sqrt(3)*0.5*v_jup, 0.5*v_jup]
m = [m_sun, m_jup, 0.01*m_jup]


#calculate initial accelerations
for j in range(0, nbodies):
    for k in range(0, nbodies):
        if j != k:
                dx = x[0, j*nbodies:j*nbodies+3] - x[0, k*nbodies:k*nbodies+3]
                r = np.sqrt(sum(dx*dx))
                a[0, j*nbodies:j*nbodies+3] = a[0, j*nbodies:j*nbodies+3] - G*m[k]*dx/abs(r**3)

#integration loop
for i in range(0, nsteps): 
    #adjust initial velocity to half timesteps
    if (i == 0):
        for j in range(0, nbodies):
            v[i, j*nbodies:j*nbodies+3] = v[0, j*nbodies:j*nbodies+3] + 0.5*dt*a[i, j*nbodies:j*nbodies+3]
    #update position
    else:
        for j in range(0, nbodies):
            x[i, j*nbodies:j*nbodies+3] = x[i-1, j*nbodies:j*nbodies+3] + dt*v[i-1, j*nbodies:j*nbodies+3]
        #calculate accelerations
        for j in range(0, nbodies):
            for k in range(0, nbodies):
                if j != k:
                    dx = x[i, j*nbodies:j*nbodies+3] - x[i, k*nbodies:k*nbodies+3]
                    r = np.sqrt(sum(dx*dx))
                    a[i, j*nbodies:j*nbodies+3] = a[i, j*nbodies:j*nbodies+3] - G*m[k]*dx/abs(r**3)
        #calculate velocity
        if (i == nsteps-1):
            v[i, j*nbodies:j*nbodies+3] = v[i-1, j*nbodies:j*nbodies+3] + 0.5*dt*a[i, j*nbodies:j*nbodies+3]
        else:
            for j in range(0, nbodies):
                v[i, j*nbodies:j*nbodies+3] = v[i-1, j*nbodies:j*nbodies+3] + dt*a[i, j*nbodies:j*nbodies+3]

fig = plt.figure()
plt.plot(x[:,3], x[:,4])
plt.plot(x[:,0], x[:,1])
plt.plot(x[:,6], x[:,7])
plt.xlabel("x (meters)")
plt.ylabel("y (meters)")
plt.show()