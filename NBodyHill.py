<<<<<<< HEAD
import astropy
import numpy as np
import matplotlib.pyplot as plt

####Mass in terms of solar mass (MS), distances in terms of astronomical units (AU), time in terms of years

#Initialize/declare a bunch of constants!
global G
global m_sun
global m_jup
G = 4*np.pi**2 #AU^3 MS ^-1 yr ^-2
m_sun = 1 #MS
m_jup = 0.00095 #MS
d_jup = 1 #max distance from sun
T = 2*np.pi*np.sqrt((d_jup**3)/(G*(m_sun+m_jup))) #calculated period of orbit 
tmax = 3*T #years
nsteps = 50000
nbodies = 3
dt = tmax/nsteps

#generate empty arrays for x, v, a + more innitialization
x, v, a = np.empty([nsteps, 3*nbodies]), np.empty([nsteps, 3*nbodies]), np.empty([nsteps, 3*nbodies])
m = np.empty([nbodies])
time = np.empty([nsteps])
Etotal = np.empty([nsteps])
dv = np.empty([nsteps])
m = [m_sun, m_jup, 0.00*m_jup] #masses of bodies
barycenter = (m_jup*d_jup)/(m_sun+m_jup)
r1 = barycenter
r2 = d_jup - barycenter
v0 = np.sqrt((m_sun)*G*r2/d_jup**2) #velocity for circular orbit of planet
v_sun = -m_jup*v0/m_sun #to give the system net zero momentum, I give the sun a small velocity
circ2 = np.sqrt(m_jup*G/(0.001*d_jup)) #circular orbit of moon

#initial values for x, v
x[0] = [-r1, 0.0, 0.0, r2, 0.0, 0.0, 1.001*d_jup - barycenter, 0.0, 0.0]
v[0] = [0.0, v_sun, 0.0, 0.0, v0, 0.0, 0.0, v0 + circ2, 0.0]
x0 = x[0, 3]
E0 = 0.5*v0**2 - G*(m_sun)/(r1+r2) #calculate energy
semimajor = -G*(m_sun + m_jup)/(2*E0) #calculate a from E
w = 2*np.pi*np.sqrt(G*(m_sun + m_jup)/(4*(np.pi**2)*(semimajor)**3)) #calculate angular frequency of orbit

#calculate initial accelerations
for j in range(0, nbodies):
    for k in range(0, nbodies):
        if j != k:
                dx = x[0, j*nbodies:j*nbodies+3] - x[0, k*nbodies:k*nbodies+3]
                r = np.sqrt(sum(dx*dx))
                a[0, j*nbodies:j*nbodies+3] = a[0, j*nbodies:j*nbodies+3] - G*m[k]*dx/abs(r**3)
#integration loop
rmin = 100
rmax = 0
for i in range(0, nsteps): 
    time[i] = i*dt
    #print(Etotal[i])
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
        #update final velocity back to integer timesteps
        if (i == nsteps-1):
            v[i, j*nbodies:j*nbodies+3] = v[i-1, j*nbodies:j*nbodies+3] + 0.5*dt*a[i, j*nbodies:j*nbodies+3]
        else:
            #normal calculate velocity
            for j in range(0, nbodies):
                v[i, j*nbodies:j*nbodies+3] = v[i-1, j*nbodies:j*nbodies+3] + dt*a[i, j*nbodies:j*nbodies+3]

#create new array for x and y in hill coordinates
xprime = np.empty([nsteps, 2*nbodies])

#rotate the matrix by an angle of w*t to change coordinate system to Hill
for i in range(0, len(x)):
    t = dt*i
    xprime[i] = np.array([np.cos(w*t)*x[i, 0] + np.sin(w*t)*x[i, 1], -1*np.sin(w*t)*x[i, 0] + np.cos(w*t)*x[i, 1], np.cos(w*t)*x[i, 3] + np.sin(w*t)*x[i, 4], -1*np.sin(w*t)*x[i, 3] + np.cos(w*t)*x[i, 4], np.cos(w*t)*x[i, 6] + np.sin(w*t)*x[i, 7], -1*np.sin(w*t)*x[i, 6] + np.cos(w*t)*x[i, 7]])

#plot in Newtonian
fig = plt.figure()
plt.plot(x[:,3], x[:,4])
plt.plot(x[:,0], x[:,1])
plt.title("Sun/Planet (Newtonian Coordinates)")
plt.axis('equal')
plt.xlabel("x (meters)")
plt.ylabel("y (meters)")
plt.show()

#plot in Hill
fig3 = plt.figure()
plt.plot(xprime[:,0], xprime[:,1])
plt.plot(xprime[:,2], xprime[:,3])
plt.title("Sun/Planet (Hill Coordinates)")
plt.axis('equal')
plt.xlabel("x' (meters)")
plt.ylabel("y' (meters)")
plt.show()
=======
import astropy
import numpy as np
import matplotlib.pyplot as plt
#Mass in terms of solar mass (MS), distances in terms of astronomical units (AU), time in terms of years
global G
global m_sun
global m_jup
G = 4*np.pi**2 #AU^3 MS ^-1 yr ^-2
m_sun = 1 #MS
m_jup = 0.0009543 #MS
d_jup = 5.2 #max distance from sun

T = 2*np.pi*np.sqrt((d_jup**3)/(G*(m_sun+m_jup))) #calculated period of orbit 
tmax = 10*T #years

nsteps = 500000
nbodies = 3
dt = tmax/nsteps

#generate empty arrays for x, y, z
x = np.empty([nsteps, 3*nbodies])
v = np.empty([nsteps, 3*nbodies])
a = np.empty([nsteps, 3*nbodies])
m = np.empty([nbodies])

m = [m_sun, m_jup, 0.000*m_jup]

circ = np.sqrt((m_sun+m_jup)*G/d_jup)
v_sun = -m_jup*circ/m_sun
circ2 = np.sqrt(m_jup*G/(0.001*d_jup))
#initial values for x, v
x[0] = [0.0, 0.0, 0.0, d_jup, 0.0, 0.0, 1.001*d_jup, 0.0, 0.0]
v[0] = [0.0, v_sun, 0.0, 0.0, circ, 0.0, 0.0, circ+circ2, 0.0]
m = [m_sun, m_jup, 0.000*m_jup]

E0 = 0.5*circ**2 - G*(m_sun+m_jup)/d_jup
semimajor = -G*(m_sun+m_jup)/(2*E0)
w = 2*np.pi*np.sqrt((m_sun+m_jup)/(semimajor**3))

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
        #update final velocity back to integer timesteps
        if (i == nsteps-1):
            v[i, j*nbodies:j*nbodies+3] = v[i-1, j*nbodies:j*nbodies+3] + 0.5*dt*a[i, j*nbodies:j*nbodies+3]
        else:
            #normal calculate velocity
            for j in range(0, nbodies):
                v[i, j*nbodies:j*nbodies+3] = v[i-1, j*nbodies:j*nbodies+3] + dt*a[i, j*nbodies:j*nbodies+3]

#create new array for x and y in hill coordinates
xprime = np.empty([nsteps, 2*nbodies])
#rotate the matrix by an angle of w*t to change coordinate system to Hill
for i in range(0, len(x)):
    t = dt*i
    xprime[i] = np.array([np.cos(w*t)*x[i, 0] + np.sin(w*t)*x[i, 1], -1*np.sin(w*t)*x[i, 0] + np.cos(w*t)*x[i, 1], np.cos(w*t)*x[i, 3] + np.sin(w*t)*x[i, 4], -1*np.sin(w*t)*x[i, 3] + np.cos(w*t)*x[i, 4], np.cos(w*t)*x[i, 6] + np.sin(w*t)*x[i, 7], -1*np.sin(w*t)*x[i, 6] + np.cos(w*t)*x[i, 7]])

#plot in Newtonian
fig = plt.figure()
plt.plot(x[:,3], x[:,4])
plt.plot(x[:,0], x[:,1])
#plt.plot(x[:,6], x[:,7])
plt.title("Sun/Planet (Newtonian Coordinates)")
plt.axis('equal')
plt.xlabel("x (meters)")
plt.ylabel("y (meters)")
plt.show()

#plot in Hill
fig3 = plt.figure()
plt.plot(xprime[:,0], xprime[:,1])
plt.plot(xprime[:,2], xprime[:,3])
#plt.plot(xprime[:,4], xprime[:,5])
plt.title("Sun/Planet (Hill Coordinates)")
plt.axis('equal')
plt.xlabel("x (meters)")
plt.ylabel("y (meters)")
plt.show()
>>>>>>> a7e629a14e292dba4f96ee8fda6060b661684260
