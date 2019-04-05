import kplr
from astropy.io import fits as pyfits
import scipy
from scipy import signal
import astropy
import numpy as np
import statistics as stat
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import time as ptime
import operator
import pandas as pd
np.seterr(divide='ignore', invalid='ignore')

def getkois(n, final, file):  
    print("Getting KOIs...")
    client = kplr.API()
    # Gets false positive KOIs from csv where KOIs are in the second column
    IDs = np.genfromtxt(file, dtype=None, usecols = 1, skip_header = 12, unpack= True, delimiter = ',').astype(str)
    IDs = IDs.tolist()
    kois = np.zeros(n, object)
    for i in range(n): 
        koi = float(IDs[i + final][1:])
        kois[i] = client.koi(koi)
        print("Reading %d..." % i)
    periods = np.zeros(n)
    lcss = np.zeros(n, object)
    durations = np.zeros(n)
    centroids = np.zeros(n)
    stars = np.zeros(n, object)
    for i in range(n):
        periods[i] = kois[i].koi_period
        lcss[i] = np.array(kois[i].get_light_curves(), dtype=object)
        durations[i] = kois[i].koi_duration
        centroids[i] = kois[i].koi_time0bk
        stars[i] = np.array(kois[i].star, dtype=object)
    return kois, periods, lcss, durations, centroids, stars

def formatlightcurves(koi, lcs, period, centroid):
    print("Formatting light curve...")
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
    lightcurve = pd.DataFrame({"time": time, "flux" : flux})
    lightcurve = lightcurve.dropna()
    centroid = np.mod(centroid, period)
    angular_frequency = (2*np.pi)/period
    # Scales the period to 2pi
    lightcurve = lightcurve.multiply(angular_frequency, axis = 'time')
    centroid = centroid*angular_frequency
    # Shifts the times so that the main transit is at pi
    shift = np.pi - centroid
    if shift>0:
        shifting = lightcurve['time'].values
        shifting = [i + shift for i in shifting]
        shifting = [j - 2*np.pi if j > 2*np.pi else j for j in shifting]
        lightcurve['time'] = shifting
    elif shift<0:
        shifting = lightcurve['time'].values
        shifting = [i + shift for i in shifting]
        shifting = [j + 2*np.pi if j < 0 else j for j in shifting]
        lightcurve['time'] = shifting

    lightcurve = lightcurve.sort_values(by=['time'])
    time = lightcurve['time'].values
    flux = lightcurve['flux'].values
    flux = flux/np.median(flux)
    return(flux, time)

def makeplots(flux, time, star):
    fig, ax = plt.subplots()
    ax.scatter(time, flux, alpha = 0.2, s = 0.8)
    ax.vlines([np.pi, 2*np.pi/3, 4*np.pi/3], min(flux), max(flux))
    ax.set_title(str(koi.star))
    ax.set_xlabel("Phase Angle")
    ax.set_ylabel("Relative Flux")
    ax.set_xlim([min(time), max(time)])
    ax.set_ylim([min(flux), max(flux)])
    return fig, ax

def filter(flux, time, duration, period):
    print("Filtering light curve...")
    duration = duration/24
    angular_frequency = 2*np.pi/period
    centroid = np.pi
    duration = duration*angular_frequency
    transit_start = centroid - duration/2
    transit_end = centroid + duration/2
    for i in range(len(time)-1):
        if (time[i] <= transit_start) and (time[i+1] >=transit_start):
            start = i
        elif (time[i] <= transit_end) and (time[i+1] >=transit_end):
            end = i
    template_flux = flux[start:end]
    template_time = time[start:end]
    fir = template_flux[::-1]
    template_length = int(len(fir)/2)
    det = signal.correlate(flux, fir)
    det_length = len(det)
    #det = signal.lfilter(fir, 1, flux)
    det = det[template_length:template_length + len(flux)]
    det = det/np.median(det)
    return det

def detect(det, time):
    print("Detecting transit...")
    curve_filtered = pd.DataFrame({"time" : time, "det" : det})
    count = 0
    min_width = int(len(det)*0.03)
    print(min_width)
    flux_minima = signal.argrelextrema(det, np.less, order = min_width)
    minima = curve_filtered.ix[flux_minima[0]]
    print(minima)
    #minima = curve_filtered[curve_filtered['det'].isin(flux_minima[0])]
    minima_times = minima['time'].tolist()
    return minima_times

if __name__ == '__main__':
    n = 8
    file = "cumulative.csv"
    minima = np.array([])
    for j in range(6):
        final = j*n
        kois, periods, lcss, durations, centroids, stars = getkois(n, final, file)
        progress = 0
        for i in range(n):
            koi = kois[i]
            lcs = lcss[i]
            period = periods[i]
            duration = durations[i]
            centroid = centroids[i]
            star = stars[i]
            flux, time = formatlightcurves(koi, lcs, period, centroid)
            fig, ax = makeplots(flux, time, star)
            plt.show()
            plt.close()
            det = filter(flux, time, duration, period)
            fig, ax = makeplots(det, time, star)
            plt.show()
            plt.close()
            curve_min = np.array(detect(det, time))
            minima = np.concatenate((minima, curve_min))    
            print(minima)      
            progress += 1
            complete = 100*progress/n
            print ("%.2f%% complete" % complete)
        nbins = int(len(minima)/2)
        #plt.hist(minima, bins = nbins)
        #plt.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
        #plt.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
        #plt.show()
        #plt.close()
    plt.hist(minima, bins = nbins)
    plt.show()
    plt.close()