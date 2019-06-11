from Filtering import *
import kplr
from astropy.io import fits as pyfits
import scipy
import astropy
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def findparabola(x, y):
	x1, x2, x3 = x
	y1, y2, y3 = y
	
	A = ((x1-x2)*(y1-y3)-(x1-x3)*(y1-y2))/((x1**2-x3**2)*(x1-x2)-(x1**2-x2**2)*(x1-x3))
	B = ((y1-y2)-A*(x1**2-x2**2))/(x1-x2)
	C = y1 - A*x1**2 - B*x1
	
	return A, B, C

def find3maxes(l):
	y1 = max(l)
	x1 = locations[np.argmax(l)]

	y2 = max(n for n in l if n != y1)
	x2 = locations[np.where(l == y2)[0]]

	y3 = max(n for n in l if n != y2 and n != y1)
	x3 = locations[np.where(l == y3)[0]]

	return np.array([x1, x2, x3]), np.array([y1, y2, y3])

curveno = 0

lightcurves = getkois('cumulative.csv', 4, 0)
lightcurves[curveno].format()
my_curve = lightcurves[curveno]
#my_curve.plot()
widths = np.array([1/4, 1/8, 1/16, 1/32]) #magic
epsilon = 0.2 #magic
maxCC = 0
rotated = my_curve.rotate()
#rotated.plot()

#find best location/width combo from grid
done = False
fig, ax = plt.subplots(figsize=(8, 6))
for w in widths:
	if done == False:
		print(f'width = {w}')
		dphi = w/5 #magic
		locations = np.arange(min(rotated.time) + epsilon, max(rotated.time) -epsilon, dphi)
		cc_list = np.zeros(len(locations))
		i = 0
		for phi in locations:
			print(f'phi is {phi}')
			my_filter = Filter(rotated, w, phi)
			my_filter.format()
			cc = my_filter.apply(rotated.flux)
			cc_list[i] = cc
			if cc>maxCC:
				maxCC = cc
				phimax = phi
				wmax = w
			i += 1
		if max(cc_list)<maxCC:
			done = True

print(f"max response width is {wmax}, location is {phimax}, response is {maxCC}")

#find ideal location for best width
i = 0
dphi = wmax/5
locations = np.arange(min(rotated.time) + epsilon, max(rotated.time) -epsilon, dphi)
cclist = np.zeros(len(locations))
for phi in locations:
	print(f'phi is {phi}')
	maxfilter = Filter(rotated, wmax, phi)
	maxfilter.format()
	#my_filter.plot()
	cc = maxfilter.apply(rotated.flux)
	cclist[i] = cc
	i += 1

x, y= find3maxes(cclist)
A, B, C = findparabola(x, y)

#plot everything
maxloc = -B/(2*A)
xlist = np.linspace(maxloc-0.5, maxloc+0.5, 1000)
parab = A*(xlist**2)+B*xlist+C
maxval = A*(maxloc**2) + B*maxloc + C
print(maxloc, maxval)
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(xlist, parab, color = 'm')
ax.set_xlim(min(rotated.time), max(rotated.time))
ax.set_ylim(min(cclist), max(cclist))
ax.axvline(x=maxloc)
ax.scatter(locations, cclist)
ax.set_xlabel('Filter Location')
ax.set_ylabel('Cross Correlation')
ax.set_title(f'Best Width: {wmax}, Best Location {maxloc}')
plt.show()
#plt.savefig(str(my_curve.star)[1:-1]+ 'response')