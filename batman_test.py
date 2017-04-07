import batman
import numpy as np
import matplotlib.pyplot as plt

params = batman.TransitParams()       #object to store transit parameters
params.t0 = 0.                        #time of inferior conjunction
params.per = 30.360446671                       #orbital period
params.rp = 0.0989                       #planet radius (in units of stellar radii)
params.a = 37.093                        #semi-major axis (in units of stellar radii)
params.inc = 89.5                      #orbital inclination (in degrees)
params.ecc = 0.                       #eccentricity
params.w = 90.                        #longitude of periastron (in degrees)
params.limb_dark = "nonlinear"        #limb darkening model
params.u = [0.5, 0.1, 0.1, -0.1]      #limb darkening coefficients

t = np.linspace(-3, 3, 10000)  #times at which to calculate light curve
m = batman.TransitModel(params, t)    #initializes model
flux = m.light_curve(params)                    #calculates light curve

radii = np.linspace(0.09, 0.11, 20)
for r in radii:
        params.rp = r                           #updates planet radius
        new_flux = m.light_curve(params)        #recalculates light curve

plt.scatter(t, flux, alpha = 0.3, s = 0.8)
plt.xlabel("Time from central transit (days)")
plt.ylabel("Relative flux")
plt.ylim((0.989, 1.001))
plt.savefig("lc.png")
plt.show()
