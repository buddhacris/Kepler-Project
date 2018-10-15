import astropy
import numpy as np
import matplotlib.pyplot as plt
import time

global G
G = 4*np.pi**2 #AU^3 MS ^-1 yr ^-2

def initialize(tmax, nsteps):  #Initialize/declare a bunch of constants!
	m1 = 1 #MS
	m2 = 1e-4 #MS
	d = 1 #max distance from sun
	T = 2*np.pi*np.sqrt((d**3)/(G*(m1+m2))) #calculated period of orbit 
	tmax = tmax*T #years
	nsteps = nsteps
	nbodies = 4
	dt = tmax/nsteps
	e = 0
	x, v = np.empty([nsteps, 3*nbodies]), np.empty([nsteps, 3*nbodies])
	m = np.empty([nbodies])
	t = np.empty([nsteps])
	m = [m1, m2, 1e-6*m2, 1e-6*m2] #masses of bodies
	r1 = (m2*d)/(m1+m2) #calculate semi-major axis from barycenter of two main bodies
	r2 = d - r1 #semi-minor axis
	v0 = np.sqrt(1-e)*np.sqrt((m1)*G*r2/d**2) #velocity for circular orbit of planet
	v_sun = -m2*v0/m1 #to give the system net zero momentum, I give the sun a small velocity
	x[0] = [-r1, 0.0, 0.0, r2, 0.0, 0.0, 0.5*d - r1, 0.5*np.sqrt(3)*d, 0.0, 0.5*d - r1, -0.5*np.sqrt(3)*d, 0.0]
	v[0] = [0.0, v_sun, 0.0, 0.0, v0, 0.0, -np.sqrt(3)*v0/2, v0/2, 0.0, np.sqrt(3)*v0/2, v0/2, 0.0]
	x0 = x[0, 3]
	E0 = 0.5*v0**2 - G*(m1)/(r1+r2) #calculate energy
	semimajor = -G*(m1 + m2)/(2*E0) #calculate a from E
	w = 2*np.pi*np.sqrt(G*(m1 + m2)/(4*(np.pi**2)*(semimajor)**3)) #calculate angular frequency of orbit

	return tmax, nsteps, nbodies, dt, e, x, v, m, t, w

def perturbation(x, v, coord, k):
	v0 = v[0, 4]
	if coord == "x":
		p = 6
		x[0, p] = x[0, p] + k*v0
		x[0, p+3] = x[0, p+3] + k*v0
	elif coord == "y":
		p = 7
		x[0, p] = x[0, p] + k*v0
		x[0, p+3] = x[0, p+3] + k*v0
	elif coord == "z":
		p = 8
		x[0, p] = x[0, p] + k*v0
		x[0, p+3] = x[0, p+3] + k*v0
	if coord == "vx":
		p = 6
		v[0, p] = v[0, p] - perturbation*v0
		v[0, p+3] = v[0, p+3] + perturbation*v0
	elif coord == "vy":
		p = 7
		v[0, p] = v[0, p] - perturbation*v0
		v[0, p+3] = v[0, p+3] + perturbation*v0
	elif coord == "vz":
		p = 8
		v[0, p] = v[0, p] - perturbation*v0
		v[0, p+3] = v[0, p+3] + perturbation*v0

	return x, v
        
def leapfrog(x, v, m, tmax, nsteps, nbodies, dt, t):
    #integration loop
    for i in range(0, nsteps): 
        t[i] = i*dt
        acceleration = [0, 0, 0, 0]  #calculate initial accelerations
        for j in range(0, nbodies):
            for k in range(0, nbodies):
                if j != k:
                    dx = x[0, j*3:j*3 + 3] - x[0, k*3:k*3 + 3]
                    r = np.sqrt(sum(dx*dx))
                    acceleration[j] = -G*m[k]*dx/abs(r**3)
        #adjust initial velocity to half timesteps
        if (i == 0):
            for j in range(0, nbodies):
                v[i, j*3:j*3 + 3] = v[0, j*3:j*3 + 3] + 0.5*dt*acceleration[j]
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
                        acceleration[j] = acceleration[j] - G*m[k]*dx/abs(r**3)
            #update final velocity back to integer timesteps
            if (i == nsteps-1):
                v[i, j*3:j*3 + 3] = v[i-1, j*3:j*3 + 3] + 0.5*dt*acceleration[j]
            else:
                #normal calculate velocity
                for j in range(0, nbodies):
                    v[i, j*3:j*3 + 3] = v[i-1, j*3:j*3 + 3] + dt*acceleration[j]

    return x

def rotate(x, w, dt):
	nbodies = int(len(x[0, :])/3)
	for i in range(0, nbodies):
		x[:, i] = x[:,i] - x[:, i%3] #finds distance from each body to sun

	#rotate the matrix by an angle of w*t to change coordinate system to Hill
	for i in range(len(x)):
		state = x[i]
		newstate = np.zeros(len(state))
		for j in range(3*nbodies):
			t = dt*i
			if (j%3 == 0):
				newstate[j] = np.cos(w*t)*state[j] + np.sin(w*t)*state[j+1]
			if (j%3 == 1):
				newstate[j] = -1*np.sin(w*t)*state[j-1] + np.cos(w*t)*state[j]
			if (j%3 == 2):
				newstate[j] = state[j]
		x[i] = newstate
	return x

def plot(x, k, coord):
	plt.plot(x[:,0], x[:,1], "y.", markersize = 1)
	plt.plot(x[:,3], x[:,4], markersize = 12)
	plt.plot(x[:,6], x[:,7], markersize = 12)
	plt.plot(x[:, 9], x[:, 10], markersize = 12)
	plt.title("Hill Coordinates: %.2E%% %s perturbation" % (k, coord))
	plt.axis('equal')
	plt.xlabel("x' (AU)")
	plt.ylabel("y' (AU)")
	#plt.show()
	plt.savefig("%sPerturbTest%.2E.png" % (coord, k))

def simulate(coord, k, scale):
	tmax, nsteps, nbodies, dt, e, x, v, m, t, w = initialize(15*scale, 40000*scale)
	x, v = perturbation(x, v, coord, k)
	x = leapfrog(x, v, m, tmax, nsteps, nbodies, dt, t)
	xprime = rotate(x, w, dt)
	plot(xprime, k, coord)

if __name__ == '__main__':
	k = 1e-4
	scales = [100]
	for scale in scales:
		start_time = time.process_time()
		simulate("x", k, scale)
		tim = (time.process_time() - start_time)
		print(tim)