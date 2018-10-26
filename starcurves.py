import kplr
import pyfits
import astropy
import numpy as np
import statistics as stat
import matplotlib.pyplot as plt
np.seterr(divide='ignore', invalid='ignore')

def getkois(file):  
    # Gets false positive KOIs from csv where KOIs are in the second column
    kois = np.genfromtxt(file, dtype=None, usecols = 1, skip_header = 12, unpack= True, delimiter = ',').astype(str)
    kois = kois.tolist()
    for i in range(len(kois)): 
        #convert string to float
        kois[i] = float(kois[i][1:])
    return kois

def formatlightcurves(koi):
    # Get a list of light curve datasets.
    lcs = koi.get_light_curves(short_cadence=False)
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
    centroid = np.mod(koi.koi_time0bk, period)
    angular_frequency = (2*np.pi)/period
    # Scales the period to 2pi
    for t in range(len(time)):
            time[t] = (time[t])*angular_frequency
    centroid = centroid*angular_frequency
    # Shifts the times so that the main transit is at pi
    shift = np.pi - centroid
    if shift>0:
        for t in range(len(time)):
            time[t] = time[t] + shift
            if time[t] > 2*np.pi:
                time[t] = time[t] - 2*np.pi
    elif shift<0:
        for t in range(len(time)):
            time[t] = time[t] + shift
            if time[t] <= 0:
                time[t] = time[t] + 2*np.pi
    return(flux, time)

def makeplots(flux, time, koi):
    fig, ax = plt.subplots()
    ax.scatter(time, flux, alpha = 0.2, s = 0.8)
    ax.vlines([np.pi, 2*np.pi/3, 4*np.pi/3], min(flux), max(flux))
    ax.set_title(str(koi.star))
    ax.set_xlabel("Phase Angle")
    ax.set_ylabel("Relative Flux")
    ax.set_xlim([min(time), max(time)])
    ax.set_ylim([min(flux), max(flux)])


if __name__ == '__main__':
    client = kplr.API()
    file = "cumulative.csv"
    kois = getkois(file)
    for ID in kois:
        koi = client.koi(ID)
        star = koi.star
        starkois = star.kois
        for koi in starkois:
            flux, time = formatlightcurves(koi)
            makeplots(flux, time, koi)
            plt.show()