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
# Plot!
plt.scatter(time, flux, alpha = 0.1, s = 2)
plt.title(str(ID))
plt.xlabel("Time (BJD)")
plt.ylabel("Flux")
plt.xlim([0,31])
plt.ylim([0.975,1.005])

plt.show()
