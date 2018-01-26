import astropy
import numpy as np
import matplotlib.pyplot as plt

m_jup = 1.898e27 #kg
d = 8.157e11 #max distance from sun
dmin = 7.409e11 #min distance from jupiter to sun
m_sun = 1.989e30 #kg
G = 6.674e-11 #m^3 kg ^-1 s ^-2
v_jup = 12.44e3 #ms^-1

semimajor = 0.5*(d+dmin) #length of semimajor axis in m
T = 2*np.pi*np.sqrt((semimajor**3)/(G*m_sun)) #calculated period of orbit 

tmax = T #seconds
count = 1000000 #timesteps
dt = tmax/count

#initial position of jupiter
x = [d]
y = [0.0]

u = [0.0] #velocity in the y-direction
v = [v_jup] #velocity in the x-direction

def integrate(x, y, u, v, G, M, m, dt, count):

	for i in range(0, count):
		r = np.sqrt(x[i]**2 + y[i]**2) #distance to sun

		ax = -G*M*x[i]/(r**3) #acceleration
		ay = -G*M*y[i]/(r**3)

		u.append( u[i] + dt*ax ) #update velocity
		v.append( v[i] + dt*ay )

		x.append( x[i] + dt*u[i] ) #update position
		y.append( y[i] + dt*v[i] )

	energy0 = 0.5*m*v[0]**2 - G*M*m/r
	energy = 0.5*m*v[-1]**2 - G*M*m/r
	print(energy0, energy)



	plt.plot(x, y)
	plt.xlabel("x (meters)")
	plt.ylabel("y (meters)")
	plt.show()

integrate(x, y, u, v, G, m_sun, m_jup, dt, count)
