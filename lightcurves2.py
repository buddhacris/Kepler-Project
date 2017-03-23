import kplr
import pyfits
import astropy
import numpy as np
import statistics as stat
import matplotlib.pyplot as plt

client = kplr.API()

# Find the target KOI.
ID = 189.01
koi = client.koi(ID)

star = koi.star
# Get a list of light curve datasets.
lcs = star.get_light_curves(short_cadence=False)

period = koi.koi_period

# Loop over the datasets and read in the data.
time, flux = [], []
for lc in lcs:
    with lc.open() as f:
        hdu_data = f[1].data
        newtime = hdu_data["time"].tolist()
        newflux = hdu_data["pdcsap_flux"].tolist()
        flux_median = stat.median(newflux)
# Dividing out the median from the flux to keep everything near 1
        for f in range(0, len(newflux)):
            newflux[f] = newflux[f]/flux_median
# Getting the mod of the BJD to stack every transit
        for t in range(0, len(newtime)):
            newtime[t] = np.mod(newtime[t], period)
# Combines current dataset with previous time and flux arrays
        flux = flux + newflux
        time = time + newtime

# Calculates centroid of transit
xcentroid = np.mod(koi.koi_time0bk, period)

time0 = min(time)
timez = max(time)

period = max(time)-min(time)

mtime = ((xcentroid - time0) - (timez - xcentroid))/2
#Re-scales so that the period is 2pi and the transit is at pi
for t in range(0, len(time)):
    if time[t] < mtime + time0:
        time[t] = time[t] + period

time0 = min(time)
timez = max(time)
angular_frequency = (2*np.pi)/period

for t in range(0, len(time)):
    time[t] = (time[t] - time0)*angular_frequency

period = max(time) - min(time)
xcentroid = (xcentroid-time0) *angular_frequency

# Plots time vs flux as scatter plot
plt.scatter(time, flux, alpha = 0.1, s = 0.8)
# Plots centroid of transit as vertical line, with two other lines 1/6 period ahead and behind
plt.axvline(x = xcentroid, alpha = 0.5)
plt.axvline(x = xcentroid + (np.pi/3), alpha = 0.4)
plt.axvline(x = xcentroid - (np.pi/3), alpha = 0.4)
# Names plot after the KOI ID
plt.title(str(koi.kepoi_name) + ", ID: " + str(koi.kepid) + ", Period: " + str(koi.koi_period))
plt.xlabel("Phase Angle")
plt.ylabel("Flux")
# fits plot x-range to transit period, y-range is specific to this KOI (will change soon)
plt.xlim([min(time), max(time)])
plt.ylim([0.975, 1.005])
# Made the plot a little bit wider
plt.figure(num=1, figsize=(20, 10))

plt.show()
