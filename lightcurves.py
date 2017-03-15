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
            newtime[t] = np.mod(newtime[t], koi.koi_period)
# Combines current dataset with previous time and flux arrays
        flux = flux + newflux
        time = time + newtime
# Calculates centroid of transit
xcentroid = np.nanmean(time)
# Plots time vs flux as scatter plot
plt.scatter(time, flux, alpha = 0.1, s = 0.8)
# Plots centroid of transit as vertical line, with two other lines 1/6 period ahead and behind
plt.axvline(x = xcentroid, alpha = 0.5)
plt.axvline(x = xcentroid + (koi.koi_period/6), alpha = 0.4)
plt.axvline(x = xcentroid - (koi.koi_period/6), alpha = 0.4)
# Names plot after the KOI ID
plt.title(str(ID))
plt.xlabel("Time (BJD)")
plt.ylabel("Flux")
# fits plot x-range to transit period, y-range is specific to this KOI (will change soon)
plt.xlim([0, koi.koi_period])
plt.ylim([0.975, 1.005])
# Made the plot a little bit wider
plt.figure(num=1, figsize=(20, 10))

plt.show()
