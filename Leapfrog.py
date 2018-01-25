import astropy
import numpy as np
import matplotlib.pyplot as plt

m_jup = 1.898e27 #kg
d = 8.157e11 #m
m_sun = 1.989e30 #kg
G = 6.674e-11 #m^3 kg ^-1 s ^-2
v_jup = 1.31e4 #ms^-1
dt = 10000 #seconds

#initial position of jupiter
x = [d]
y = [0.0]

u = [0.0] #velocity in the y-direction
v = [v_jup] #velocity in the x-direction

def integrate(x, y, u, v, G, M, dt):
	for i in range(0, 100000):
		r = np.sqrt(x[i]**2 + y[i]**2)

		ax = G*M*x[i]/(r**3)
		ay = G*M*y[i]/(r**3)

		u.append( u[i] - dt*ax )
		v.append( v[i] - dt*ay )

		x.append( x[i] + dt*u[i] )
		y.append( y[i] + dt*v[i] )

	plt.plot(x, y)
	plt.xlabel("x (meters)")
	plt.ylabel("y (meters)")
	plt.show()

integrate(x, y, u, v, G, m_sun, dt)
