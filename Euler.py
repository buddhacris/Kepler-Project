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


#initial values for x, v

def integrate(dmax, v_jup, tmax, nsteps):
	dt = tmax/nsteps

	#generate empty arrays for x, y, z
	x = np.empty([nsteps, 3])
	v = np.empty([nsteps, 3])


	x[0] = [dmax, 0.0, 0.0]
	v[0] = [0.0, v_jup, 0.0]

	#save initial distance and velocity
	r0 = np.sqrt(sum(x[0]*x[0]))
	v0 = np.sqrt(sum(v[0]*v[0]))

	a = np.empty([nsteps, 3])
	a[0] = -G*m_sun*x[0]/(r0**3)

	#calculate initial energy and angular momentum
	Energy0 = 0.5*m_jup*v0**2 - G*m_sun*m_jup/r0
	L0 = np.cross(x[0], m_jup*v[0])

	for i in range(1, nsteps):
		v[i] = v[i-1] + dt*a[i-1]
		x[i] = x[i-1] + dt*v[i-1]

		r = np.sqrt(sum(x[i]*x[i]))
		a[i] = -G*m_sun*x[i]/(r**3)

	#calculates final distance and velocity
	rf = np.sqrt(sum(x[-1]*x[-1]))
	vf = np.sqrt(sum(v[-1]*v[-1]))

	Energyf = 0.5*m_jup*vf**2 - G*m_sun*m_jup/rf
	Lf = np.cross(x[-1], m_jup*v[-1])

	dE = abs((Energyf - Energy0)/Energy0)
	dL = abs((Lf[2]-L0[2])/L0[2])

	return dE, dL


tmax = 8*T #seconds

sizes = [10000, 50000, 100000, 500000, 1000000]
Energies = np.empty([len(sizes)])
Momenta = np.empty([len(sizes)])


for i in range(0, len(sizes)):
	nsteps = sizes[i]
	Energies[i], Momenta[i] = integrate(dmax, v_jup, tmax, nsteps)

dts = tmax/sizes
plt.scatter(dts, Energies)
#plt.yscale('log')
#plt.xscale('log')
plt.title("Euler")
plt.xlabel("Step size (seconds)")
plt.ylabel("Fractional Energy Difference")
plt.show()

plt.scatter(dts, Momenta)
#plt.yscale('log')
#plt.xscale('log')
plt.title("Euler")
plt.xlabel("Step size (seconds)")
plt.ylabel("Fractional Angular Momentum Difference")
plt.show()