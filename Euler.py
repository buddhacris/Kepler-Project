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

tmax = 8*T #seconds
count = 1000000 #timesteps
dt1 = tmax/count

count2 = count*2 #calculates time steps half the size of the original
dt2 = tmax/count2

count3 = count*2 #calculates timesteps double the size of th original
dt3 = tmax/count3

#initial position of jupiter
x = [d]
y = [0.0]

u = [0.0] #velocity in the y-direction
v = [v_jup] #velocity in the x-direction

def integrate(x, y, u, v, G, M, m, dt, count):

	r = np.sqrt(x[0]**2 + y[0]**2) #distance to sun
	energy0 = 0.5*m*v[0]**2 - G*M*m/r  #calculates initial energy

	L0 = [m*u[0]*x[0], m*v[0]*y[0]] #calculates initial angular momentum vector

	for i in range(0, count):
		r = np.sqrt(x[i]**2 + y[i]**2) #distance to sun

		ax = -G*M*x[i]/(r**3) #acceleration
		ay = -G*M*y[i]/(r**3)

		u.append( u[i] + dt*ax ) #update velocity
		v.append( v[i] + dt*ay )

		x.append( x[i] + dt*u[i] ) #update position
		y.append( y[i] + dt*v[i] )

	energy = 0.5*m*v[count-1]**2 - G*M*m/r #calculates final energy
	Lf = [m*u[count-1]*x[count-1], m*v[count-1]*y[count-1]] #calculates final angular momentum vector
	dL = [Lf[0]-L0[0], Lf[1]-L0[1]] #calculates the difference in angular momentum vectors
	L = np.sqrt(dL[0]**2 + dL[1]**2) #takes the magnitude of the previous difference vector

	"""
	print("Initial position (m): (", x[0], ",", y[0], "); Final position: (", x[-1], ",", y[-1], ")")
	print("Initial velocity (m): (", u[0], ", ", v[0], "); Final velocity: (", u[-1], ",", v[-1], ")")	
	print("Step size (sec): ", dt, "; number of steps: ", count)
	print("Initial energy (J): ", energy0, "; Final energy: ", energy)
	"""



	"""
	plt.plot(x, y)
	plt.xlabel("x (meters)")
	plt.ylabel("y (meters)")
	plt.show()
	"""
	dE = energy-energy0
	print("Energy difference: ", dE, "\nAngular momentum difference: ", L)
	return dE, L

print("Half dt:\n")
dE1, L1 = integrate(x, y, u, v, G, m_sun, m_jup, dt2, count2)
print("Original dt:\n")
dE2, L2 = integrate(x, y, u, v, G, m_sun, m_jup, dt1, count)
print("Double dt:\n")
dE3, L3 = integrate(x, y, u, v, G, m_sun, m_jup, dt3, count3)


plt.plot([dt1, dt2, dt3], [dE1, dE2, dE3])
plt.title("Euler")
plt.xlabel("Step size (seconds)")
plt.ylabel("Energy Difference")
plt.show()

plt.plot([dt1, dt2, dt3], [L1, L2, L3])
plt.title("Euler")
plt.xlabel("Step size (seconds)")
plt.ylabel("Angular Momentum Difference")
plt.show()