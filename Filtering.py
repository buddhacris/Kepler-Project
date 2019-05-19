import kplr
from astropy.io import fits as pyfits
import scipy
import astropy
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def find_nearest(array, value):
    ''' Find nearest value in an array '''
    idx = (np.abs(array-value)).argmin()
    return idx

class Curve:
    
    def __init__(self, period, duration, centroid, star, lc_list):
        self.period = period
        self.duration = duration
        self.centroid = centroid
        self.star = star
        self.lc_list = lc_list
        
        #finds time and flux upon initialization
        time, flux = [], []
        for lc in self.lc_list:
            with lc.open() as f:
                hdu_data = f[1].data
                time = time + hdu_data["time"].tolist()
                flux = flux + hdu_data["pdcsap_flux"].tolist()
            lightcurve = pd.DataFrame({"time": time, "flux" : flux})
            #drop all nans from the lightcurve
            lightcurve = lightcurve.dropna()
        self.time = lightcurve['time'].values
        self.flux = lightcurve['flux'].values
        
    def set_time(self, time):
        self.time = time
    def set_flux(self, flux):
        self.flux = flux
    def set_period(self, period):
        self.period = period
    def set_duration(self, duration):
        self.duration = duration
    def set_centroid(self, centroid):
        self.centroid = centroid
        
    def format(self):
        if hasattr(self, 'flux'):
            #this will normalize the flux of a lightcurve, wrap the time over the period, and put the centroid of the primary transit at pi
            print("Formatting light curve...")
            lightcurve = pd.DataFrame({'time': self.time, 'flux': self.flux})
            #normalize the lightcurve by dividing out the median
            med = lightcurve['flux'].median()
            prd = self.period
            lightcurve, angular_frequency = fold(lightcurve, prd, med)
            #shifts the lightcurve such that the primary transit centroid is at pi
            newcentroid = np.mod(self.centroid, prd)*angular_frequency
            lightcurve = shift(lightcurve, newcentroid)
            #sorts lightcurve by time for easy operations later
            lightcurve = lightcurve.sort_values(by=['time'])
            #updates original curve with new normalized flux, 
            self.time = lightcurve['time'].values
            self.flux = lightcurve['flux'].values
            self.centroid = newcentroid
        else:
            raise AttributeError('Lightcurve has no attribute \'flux\'')
    
    def lplot(self):
        #plots lightcurve and vertical lines at lagrange points
        fig, ax = plt.subplots()
        ax.scatter(self.time, self.flux, alpha = 0.2, s = 0.8)
        ax.vlines([np.pi, 2*np.pi/3, 4*np.pi/3], min(self.flux), max(self.flux))
        ax.set_title(str(self.star))
        ax.set_xlabel("Phase Angle [radians]")
        ax.set_ylabel("Normalized Flux [$e^-$ $s^{-1}$]")
        ax.set_xlim([min(self.time), max(self.time)])
        ax.set_ylim([min(self.flux), max(self.flux)])
        plt.show()
    
    def plot(self):
        #plots lightcurve
        fig, ax = plt.subplots()
        ax.scatter(self.time, self.flux, alpha = 0.2, s = 0.8)
        ax.set_title(str(self.star))
        ax.set_xlabel("Phase Angle [radians]")
        ax.set_ylabel("Normalized Flux [$e^-$ $s^{-1}$]")
        ax.set_xlim([min(self.time), max(self.time)])
        ax.set_ylim([min(self.flux), max(self.flux)])
        plt.show()

class Filter:
    def __init__(self, curve, width = None, location = None):
        #creates filter based off of lightcurve, with specified width and location if included
        t = curve.time
        f = curve.flux
        angular_frequency = (2*np.pi)/curve.period
        transit_duration = angular_frequency*(curve.duration/24)
        transit_start = curve.centroid - transit_duration/2
        transit_end = curve.centroid + transit_duration/2
        t1 = find_nearest(t, transit_start)
        t2 = find_nearest(t, transit_end)
        if width is None or location is None:
            filter_centroid = np.pi
            filter_width = transit_duration
            self.width = filter_width
            self.location = filter_centroid
        else:
            self.width = width
            self.location = location
        flux = np.full(3*len(t), transit_duration/(2*np.pi))
        for i in range(len(t)):
            if (t[i] > self.location - self.width/2) and (t[i] < self.location + self.width/2):
                flux[i + len(t)] = self.width/(2*np.pi)-1
        self.flux = flux
        self.time = t
        self.curve = curve
        
    def plot(self):
        plt.title("Box Matched Filter")
        plt.xlabel("Phase angle [radians]")
        plt.ylabel("Normalized Flux [$e^-$ $s^{-1}$]")
        plt.plot(self.time, self.flux[len(self.time):2*len(self.time)])
        plt.show()
        
    def apply(self):
        #applies filter to model curve
        DV = 3
        hop_size = self.width/DV
        n_hops = int(2*np.pi//hop_size)
        hop_slice_width = len(self.curve.flux)//(n_hops)
        icentroid = find_nearest(self.curve.time, np.pi)
        filter_centroid = len(self.curve.time) + icentroid
        product = np.zeros(n_hops)
        dot_times = np.zeros(n_hops)
        for i in range(n_hops):
            j = i*hop_slice_width
            slice = self.flux[(filter_centroid-j):(filter_centroid-j+len(self.curve.time))]
            dot_times[i] = self.curve.time[j]
            product[i] = np.dot(slice, self.curve.flux)
        product_centroid = find_nearest(dot_times, np.pi)
        result = product[product_centroid]
        return dot_times, product

def getkois(file, groupsize, igroup):
    print("Getting KOIs...")
    client = kplr.API()
    # Gets false positive KOIs from csv where KOIs are in the second column
    IDs = np.genfromtxt(file, dtype=str, usecols = 1, skip_header = 12, unpack= True, delimiter = ',').tolist()
    kois = np.zeros(groupsize, object)
    for i in range(groupsize): 
        koi = float(IDs[i + igroup][1:])
        kois[i] = client.koi(koi)
        print(f"Reading {i}...")
    lightcurves = np.empty(groupsize, dtype = Curve)
    for i in range(groupsize):
        lightcurves[i] = Curve(kois[i].koi_period, kois[i].koi_duration, kois[i].koi_time0bk, kois[i].star, np.array(kois[i].get_light_curves(), dtype=object))
    return lightcurves

def fold(lightcurve, prd, med):
    #folds lightcurve, and scales period to 2pi
    lightcurve['time'] = lightcurve['time'].mod(prd)
    lightcurve['flux'] = lightcurve['flux'].apply(lambda x: x/med)
    angular_frequency = (2*np.pi)/prd
    lightcurve['time'] = lightcurve.multiply(angular_frequency, axis = 'time')
    return lightcurve, angular_frequency

def shift(lightcurve, newcentroid):
    shift = np.pi - newcentroid #positive if centroid is left of pi, negative if centroid is right of pi
    newcentroid = newcentroid + shift
    if shift>0:
        lightcurve['time'] = lightcurve['time'].apply(lambda x: x + shift)
        lightcurve['time'] = lightcurve['time'].apply(lambda x: x - 2*np.pi if x>2*np.pi else x)
    elif shift<0:
        lightcurve['time'] = lightcurve['time'].apply(lambda x: x + shift)
        lightcurve['time'] = lightcurve['time'].apply(lambda x: x + 2*np.pi if x<0 else x)
    return lightcurve