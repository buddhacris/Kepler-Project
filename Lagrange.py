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
m_jup = 1e-4 #MS
d_jup = 1 #max distance from sun
T = 2*np.pi*np.sqrt((d_jup**3)/(G*(m_sun+m_jup))) #calculated period of orbit 
tmax = 5*T #years
nsteps = 8000
nbodies = 4
dt = tmax/nsteps
eccentricity = 0

#generate empty arrays for x, v, a + more initialization
x, v, a = np.empty([nsteps, 3*nbodies]), np.empty([nsteps, 3*nbodies]), np.empty([nsteps, 3*nbodies])
m = np.empty([nbodies])
time = np.empty([nsteps])
Etotal = np.empty([nsteps])
dv = np.empty([nsteps])
m = [m_sun, m_jup, 1e-6*m_jup, 1e-6*m_jup] #masses of bodies
barycenter = (m_jup*d_jup)/(m_sun+m_jup)
r1 = barycenter
r2 = d_jup - barycenter
v0 = np.sqrt(1-eccentricity)*np.sqrt((m_sun)*G*r2/d_jup**2) #velocity for circular orbit of planet
v_sun = -m_jup*v0/m_sun #to give the system net zero momentum, I give the sun a small velocity
circ2 = np.sqrt(m_jup*G/(0.001*d_jup)) #circular orbit of moon

#initial values for x, v
x[0] = [-r1, 0.0, 0.0, r2, 0.0, 0.0, 0.5*d_jup - r1, 0.5*np.sqrt(3)*d_jup, 0.0, 0.5*d_jup - r1, -0.5*np.sqrt(3)*d_jup, 0.0]
v[0] = [0.0, v_sun, 0.0, 0.0, v0, 0.0, -np.sqrt(3)*v0/2, v0/2, 0.0, np.sqrt(3)*v0/2, v0/2, 0.0]
x0 = x[0, 3]
E0 = 0.5*v0**2 - G*(m_sun)/(r1+r2) #calculate energy
semimajor = -G*(m_sun + m_jup)/(2*E0) #calculate a from E
w = 2*np.pi*np.sqrt(G*(m_sun + m_jup)/(4*(np.pi**2)*(semimajor)**3)) #calculate angular frequency of orbit

#calculate initial accelerations
for j in range(0, nbodies):
    for k in range(0, nbodies):
        if j != k:
            dx = x[0, j*3:j*3 + 3] - x[0, k*3:k*3 + 3]
            r = np.sqrt(sum(dx*dx))
            a[0, j*3:j*3 + 3] = a[0, j*3:j*3 + 3] - G*m[k]*dx/abs(r**3)
#integration loop
rmin = 100
rmax = 0
for i in range(0, nsteps): 
    time[i] = i*dt
    #adjust initial velocity to half timesteps
    if (i == 0):
        for j in range(0, nbodies):
            v[i, j*3:j*3 + 3] = v[0, j*3:j*3 + 3] + 0.5*dt*a[i, j*3:j*3 + 3]
    #update position
    else:
        for j in range(0, nbodies):
            x[i, j*3:j*3 + 3] = x[i-1, j*3:j*3 + 3] + dt*v[i-1, j*3:j*3 + 3]
        #calculate accelerations
        for j in range(0, nbodies):
            for k in range(0, nbodies):
                if j != k:
                    dx = x[i, j*3:j*3 + 3] - x[i, k*3:k*3 + 3]
                    r = np.sqrt(sum(dx*dx))
                    if (j == 1 and k == 0):
                        if (r > rmax):
                            rmax = r
                        if (r < rmin):
                            rmin = r
                    if (k == 1 and j == 0):
                        if (r > rmax):
                            rmax = r
                        if (r < rmin):
                            rmin = r
                    a[i, j*3:j*3 + 3] = a[i, j*3:j*3 + 3] - G*m[k]*dx/abs(r**3)
        #calculate velocity
        #update final velocity back to integer timesteps
        if (i == nsteps-1):
            v[i, j*3:j*3 + 3] = v[i-1, j*3:j*3 + 3] + 0.5*dt*a[i, j*3:j*3 + 3]
        else:
            #normal calculate velocity
            for j in range(0, nbodies):
                v[i, j*3:j*3 + 3] = v[i-1, j*3:j*3 + 3] + dt*a[i, j*3:j*3 + 3]

ecc = (rmax-rmin)/(rmax+rmin)

#create new array for x and y in hill coordinates
xprime = np.empty([nsteps, 2*nbodies])
xdiff = np.empty([nsteps, 3*nbodies])
for i in range(0, 3*nbodies):
    xdiff[:, i] = x[:,i] - x[:, i%3] #finds distance from each body to sun

#rotate the matrix by an angle of w*t to change coordinate system to Hill
for i in range(0, len(x)):
    t = dt*i
    xprime[i] = np.array([np.cos(w*t)*xdiff[i, 0] + np.sin(w*t)*xdiff[i, 1], -1*np.sin(w*t)*xdiff[i, 0] + np.cos(w*t)*xdiff[i, 1], np.cos(w*t)*xdiff[i, 3] + np.sin(w*t)*xdiff[i, 4], -1*np.sin(w*t)*xdiff[i, 3] + np.cos(w*t)*xdiff[i, 4], np.cos(w*t)*xdiff[i, 6] + np.sin(w*t)*xdiff[i, 7], -1*np.sin(w*t)*xdiff[i, 6] + np.cos(w*t)*xdiff[i, 7], np.cos(w*t)*xdiff[i, 9] + np.sin(w*t)*xdiff[i, 10], -1*np.sin(w*t)*xdiff[i, 9] + np.cos(w*t)*xdiff[i, 10]])


#plot in Newtonian
plt.subplot(1, 2, 1)
plt.plot(x[:,3], x[:,4], markersize = 12)
plt.plot(x[:,0], x[:,1], markersize = 12)
plt.plot(x[:,6], x[:,7], markersize = 12)
plt.plot(x[:, 9], x[:, 10], markersize = 12)
plt.title("Sun/Planet (Newtonian Coordinates)")
plt.axis('equal')
plt.xlabel("x (meters)")
plt.ylabel("y (meters)")

#plot in Hill
plt.subplot(1, 2, 2)
plt.plot(xprime[:,0], xprime[:,1], markersize = 12)
plt.plot(xprime[:,2], xprime[:,3], markersize = 12)
plt.plot(xprime[:,4], xprime[:,5], markersize = 12)
plt.plot(xprime[:, 6], xprime[:, 7], markersize = 12)
plt.title("Eccentricity = %f, mass ratio = %e, step size = %f years" % (eccentricity, m_jup, dt))
plt.axis('equal')
plt.xlabel("x' (meters)")
plt.ylabel("y' (meters)")
plt.show()