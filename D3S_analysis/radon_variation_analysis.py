import importlib
import io
import os
import csv
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.dates as mdates
from dateutil.parser import parse
from datetime import datetime
from datetime import timedelta

# Python 2 and 3: easiest option
from future.standard_library import install_aliases
install_aliases()
from urllib.parse import urlparse, urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError

import pytz
from matplotlib.backends.backend_pdf import PdfPages

import weather_data_tools as weather
importlib.reload(weather)
import spectra_fitting_tools as fitter
importlib.reload(fitter)

#--------------------------------------------------------------------------#
# Process input data
#--------------------------------------------------------------------------#
def make_int(lst): 
    '''
    Makes all entries of a list an integer
    '''
    y = []
    for i in lst:
        y.append(int(i))
    return y

def make_array(lst): 
    '''
    Makes list into an array. Also splices out the irrelevant stuff 
    for a spectra
    '''
    y = np.asarray(make_int(lst[12:]))
    return y

def get_times(rows, n, tstart, tstop):
    '''
    Get list of times for data: determines time as the midpoint between the upper and lower bounds in the integration window

    Arguments:
      - full list of inputs from data csv
      - number of hours to integrate for each data point
      - start/stop dates

    Returns:
      - list of times
    '''
    ndays = (tstop - tstart).days
    entries = 12*n
    nintervals = (24/n)
    i = 0
    counter = 0
    times = []
    while counter < ndays*nintervals:
        integration = rows[(i*entries)+1:((i+1)*entries)+1]
        i+=1

        time_range = []
        datatz = parse(integration[-1][1]).tzinfo
        if (parse(integration[-1][1])<tstop.replace(tzinfo=datatz)) and \
            (parse(integration[0][1])>tstart.replace(tzinfo=datatz)):
            for j in integration:
                time_range.append(parse(j[1]))
            times.append(time_range[int(len(time_range)/2)])
            counter+=1

    return times

def get_arrays(values_w_errs):
    vals = np.asarray([i[0] for i in values_w_errs])
    errs = np.asarray([i[1] for i in values_w_errs])
    return vals,errs

def varify_data(means,sigmas,amps):
    # check for bad fits and use average of surrounding good fits
    for i in range(len(means)):
        if means[i][1] > 100 or math.isnan(means[i][0]):
            print('Fit {} is bad!'.format(i))
            j = 1
            k = 1
            if i<(len(means)-j):
                while means[i+j][1] > 100:
                    j += 1
                    print('Trying {}+{} out of {}'.format(i,j,len(means)))
                    if i >= (len(means)-j):
                        print('Abort!')
                        break
            if i>k:
                while means[i-k][1] > 100 or math.isnan(means[i-k][0]):
                    k += 1
                    if i<k:
                        break
            if i>k and i<(len(means)-j):
                print('Averaging over {} and {}'.format(i-k,i+j))
                means[i][0] = (means[i+j][0]+means[i-k][0])/2.0
                means[i][1] = (means[i+j][1]+means[i-k][1])/2.0
                sigmas[i][0] = (sigmas[i+j][0]+sigmas[i-k][0])/2.0
                sigmas[i][1] = (sigmas[i+j][1]+sigmas[i-k][1])/2.0
                amps[i][0] = (amps[i+j][0]+amps[i-k][0])/2.0
                amps[i][1] = (amps[i+j][1]+amps[i-k][1])/2.0
            elif i<k and i<(len(means)-j):
                print('Using {}'.format(i+j))
                means[i][0] = means[i+j][0]
                means[i][1] = means[i+j][1]
                sigmas[i][0] = sigmas[i+j][0]
                sigmas[i][1] = sigmas[i+j][1]
                amps[i][0] = amps[i+j][0]
                amps[i][1] = amps[i+j][1]
            elif i>k and i>=(len(means)-j):
                print('Using {}'.format(i-k))
                means[i][0] = means[i-k][0]
                means[i][1] = means[i-k][1]
                sigmas[i][0] = sigmas[i-k][0]
                sigmas[i][1] = sigmas[i-k][1]
                amps[i][0] = amps[i-k][0]
                amps[i][1] = amps[i-k][1]
            else:
                print('Nothing makes sense')
    return means,sigmas,amps

def find_time_match(times,time,delta):
    first = 0
    last = len(times)-1
    found = False
    index = -1

    if not time.tzinfo:
        time = time.replace(tzinfo=times[0].tzinfo)
    while first<=last and not found:
        midpoint = int((first + last)/2)
        list_time = times[midpoint]
        if not list_time.tzinfo:
            list_time = list_time.replace(tzinfo=time.tzinfo)
        if abs(list_time-time) < delta :
            index = midpoint
            found = True
        else:
            if time < list_time:
                last = midpoint-1
            else:
                first = midpoint+1
    return index

def SelectDataTimeRange(start_time,stop_time,data,times):
    dataarray = np.array(data)
    timesarray = np.array(times)

    indices = np.where((timesarray>=start_time)&(timesarray<=stop_times))
    subdata  = dataarray[indices]
    subdatatimes = timesarray[indices]
   
    return subdata, subdatatimes

def merge_data(times1,data1,times2,data2):
    merged_data1 = []
    merged_data2 = []
    merged_times = []
    for i in range(len(times1)):
        time_index = find_time_match(times2,times1[i],timedelta(minutes=30))
        if time_index >= 0:
            merged_data1.append(data1[i])
            merged_data2.append(data2[time_index])
            merged_times.append(times1[i])
    return merged_times,merged_data1,merged_data2

def inTimeRange(time_string,tstart,tstop):
    time = tstart - timedelta(minutes=1)
    if isinstance(time_string, str):
        try:
            time = parse(time_string)
        except:
            print('{} Not a time!'.format(time_string))
            return False
    elif isinstance(time_string, datetime):
        time = time_string

    # check that tzinfo is set for tz aware comparisons
    if tstart.tzinfo==None:
        tstart = tstart.replace(tzinfo=time.tzinfo)
    if tstop.tzinfo==None:
        tstop = tstop.replace(tzinfo=time.tzinfo)
    #print('Checking {} > {} and < {} = {}'.format(time,tstart,tstop,(time > tstart and time < tstop)))
    return (time > tstart and time < tstop)

def get_spectra(rows, nhours, tstart, tstop):
    datatz = rows[-1][1].tzinfo
    date_itr = tstart
    times = []
    spectra = []
    counter = 0
    # break data up into days to speed up range selection
    while date_itr < tstop:
        next_day = date_itr+timedelta(days=1)
        daily_row = [row for row in rows if \
                        inTimeRange(row[1],date_itr,next_day)]
        time_itr = date_itr
        date_itr = next_day
        while time_itr < date_itr:
            time_next = time_itr+timedelta(hours=nhours)
            integration = [row for row in rows if \
                            inTimeRange(row[1],time_itr,time_next)]
            time_itr = time_next

            if len(integration)==0:
                continue

            array_lst = [] 
            for j in integration:
                array_lst.append(make_array(j))
            integrated = sum(array_lst)
            spectra.append(integrated)
            times.append(integration[int(len(integration)/2)][1])
    return times, spectra

def get_calibrations(spectra, fit_function, fit_args):
    counter = 0
    calibrations = []
    calibration_errs = []
    energy_spectra = []
    last_calib = 2.5 # default calibration
    last_err = 0
    for spectrum in spectra:
        mean,simga,amp = fit_function(spectrum,counter,*fit_args)
        calib = (1460)/(mean[0])
        calib_err = (1460)/(mean[0])**2*np.sqrt(mean[1]**2)
        if calib < 0 or calib > 10 or math.isnan(calib):
            print('invalid calibration {}, using {}'.format(calib,last_calib))
            calib = last_calib
            calib_err = last_err
        else:
            last_calib = calib
            last_err = calib_err
        calibrations.append(calib)
        calibration_errs.append(calib_err)
        energy_spectrum = np.array(spectrum)*calib
        energy_spectra.append(energy_spectrum)
        counter += 1
    return calibrations, calibration_errs

def calibrate_spectra(spectra, calibrations, times, nsum):
    E_spectra = []
    bin_times = []
    spectra_sum = []
    itimes = []
    isum = 0
    for i in range(len(spectra)):
        # list of energies = channel number * calibration (assume linear)
        energies = np.array(range(len(spectra[i])))*calibrations[i]
        print(energies)
        spectrum = np.zeros(600)
        for j in range(len(spectra[i])):
            count = spectra[i][j]
            # energy bin width = 5keV
            index = int(energies[j]/5)
            spectrum[index] += count
        if isum < nsum:
            spectra_sum.append(spectrum)
            itimes.append(times[i])
            isum += 1
        else:
            E_spectra.append(sum(spectra_sum))
            bin_times.append(itimes[int(len(itimes)/2)])
            itimes = []
            spectra_sum = []
            isum = 0
    return E_spectra, bin_times

def get_peak_fits(spectra, fit_function, fit_args):
    means = []
    sigmas = []
    amps = []
    counter = 0
    for spectrum in spectra:
        mean,sigma,amp = fit_function(spectrum,counter,*fit_args)
        means.append(mean)
        sigmas.append(sigma)
        amps.append(amp)
        counter += 1
    means,sigmas,amps = varify_data(means,sigmas,amps)
    return means,sigmas,amps

def get_peaks(rows, nhours, tstart, tstop, fit_function, fit_args):
    '''
    Applies double gaussian + expo fits to all data over some range of time

    Arguments:
      - full list of csv data input rows
      - number of hours to integrate each calculation over
      - start/stop times to run over
      - peak fitting method
      - arguments to be fed to the peak fitting method

    Returns:
      - lists of means,sigmas,amps from all gaussian fits
        - each entry in list includes the value and uncertainty
    '''
    datatz = rows[-1][1].tzinfo
    date_itr = tstart
    times = []
    means = []
    sigmas = []
    amps = []
    counter = 0
    # break data up into days to speed up range selection
    while date_itr < tstop:
        next_day = date_itr+timedelta(days=1)
        daily_row = [row for row in rows if \
                        inTimeRange(row[1],date_itr,next_day)]
        time_itr = date_itr
        date_itr = next_day
        while time_itr < date_itr:
            time_next = time_itr+timedelta(hours=nhours)
            integration = [row for row in rows if \
                            inTimeRange(row[1],time_itr,time_next)]
            time_itr = time_next

            if len(integration)==0:
                continue

            array_lst = [] 
            for j in integration:
                array_lst.append(make_array(j))
            integrated = sum(array_lst)
            mean,sigma,amp = fit_function(integrated,counter,*fit_args)
            counter += 1
            means.append(mean)
            sigmas.append(sigma)
            amps.append(amp)
            times.append(integration[int(len(integration)/2)][1])

    means,sigmas,amps = varify_data(means,sigmas,amps)
    return times,means,sigmas,amps

def get_weather_data(location,nhours,start_day,stop_day):
    tstart = parse(start_day)
    tstop = parse(stop_day)
    date_itr = tstart
    times = []
    temps = []
    while date_itr < tstop:
        data = weather.weather_station_data_scrape(location, date_itr)
        time_itr = date_itr
        date_itr = date_itr+timedelta(days=1)
        if not data:
            print('No weather data for {}'.format(date_itr))
        while time_itr < date_itr:
            time_next = time_itr+timedelta(hours=nhours)
            integration = [row for row in data if \
                            inTimeRange(row[0],time_itr,time_next)]
            time_itr = time_next
            if len(integration)==0:
                continue

            times.append(integration[int(len(integration)/2)][0])
            temps.append(np.mean(np.asarray([x[1] for x in integration])))

    return times,temps

def cut_outliers(array):
    mean, sigma = get_stats(array)
    for i in range(len(array)):
        if (array[i]>mean+5*sigma) or (array[i]<mean-5*sigma):
            if i > 0 and i < len(array)-1:
                array[i] = (array[i-1] + array[i+1])/2
            elif i==0:
                if (array[i+1]<mean+5*simga) and (array[i+1]>mean-5*simga):
                    array[i] = array[i+1]
                else:
                    array[i] = mean
            elif i==len(array)-1:
                array[i] = array[i-1]
    return array

def get_stats(array):
    return np.mean(array), np.sqrt(np.var(array))

def make_plot(points,data,errs,xlbl,ylbl,tstr,style,clr,ymin=0,ymax=0):
    fig, ax = plt.subplots()
    fig.patch.set_facecolor('white')
    plt.title(tstr)
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    if ymin and ymax:
        plt.ylim(ymin,ymax)
    ax.plot(points,data,style)
    ax.errorbar(points,data,yerr=errs,fmt=style,ecolor=clr)
    fig.autofmt_xdate()

def import_csv(url,start,stop):
    print(url)
    response = urlopen(url)
    reader = csv.reader(io.TextIOWrapper(response))
    rows = [row for row in reader if \
            inTimeRange(row[1],parse(start),parse(stop))]
    print('extracted {} entries from data url'.format(len(rows)))
    # remove meta data
    return rows

def select_data(rows,start_day,stop_day):
    tstart = parse(start_day)
    tstop = parse(stop_day)
    for row in rows:
        if isinstance(row[1], str):
            row[1] = parse(row[1])
    rows = [row for row in rows if \
            inTimeRange(row[1],tstart,tstop)]
    times, spectra = get_spectra(rows,1,tstart,tstop)
    return times,spectra

def main(times,spectra,nhours,stationID=0,wtimes=[],temps=[]):
    #---------------------------------------------------------------------#
    # Get fit results for ndays integrating over nhours for Potassium
    #---------------------------------------------------------------------#
    # single_peak_fit args: channel lims, expo offset, plot flag
    #args = [210,310,100,False]
    #args = [180,280,100,True]
    args = [360,780,7.0,100,False,'K']
    calibs,calib_err = get_calibrations(spectra, fitter.single_peak_fit,args)
    E_spectra, bin_times = calibrate_spectra(spectra,calibs,times,nhours)

    args = [180,390,7.0,100,False,'K']
    K_peaks, K_sigmas, K_amps = get_peak_fits(E_spectra, \
                                              fitter.single_peak_fit,args)
    #-------------------------------------------------------------------------#
    # Varify and break apart mean,sigma,amp values and uncertainties
    #-------------------------------------------------------------------------#
    K_ch, K_ch_errs = get_arrays(K_peaks)
    K_sig = [i[0] for i in K_sigmas]
    K_A = [i[0] for i in K_amps]

    K_ch_ave, K_ch_var = get_stats(K_ch)
    K_counts = fitter.get_peak_counts(K_ch,K_sig,K_A)

    K_count = cut_outliers(K_counts)
    K_mean, K_var = get_stats(np.asarray(K_counts))

    for i in range(len(K_ch)):
        if abs(K_ch[i]-K_ch_ave) > 3*K_ch_var:
            print('Bad K-40 fit: peak channel = {}'.format(K_ch[i]))

    #---------------------------------------------------------------------#
    # Do the same for Bizmuth-214
    #---------------------------------------------------------------------#
    # double_peak_fit args: channel lims, gaus index, expo offset, plot flag
    #args = [50,130,1,1,True]

    if stationID==0:
        args = [50,130,1,1,False,'Bi']
        Bi_peaks,Bi_sigmas,Bi_amps = get_peak_fits(E_spectra, \
                                                   fitter.double_peak_fit,args)
    if stationID==1:
        args = [90,150,5.0,1,False,'Bi']
        Bi_peaks,Bi_sigmas,Bi_amps = get_peak_fits(E_spectra, \
                                                   fitter.single_peak_fit,args)
    Bi_ch, Bi_ch_errs = get_arrays(Bi_peaks)
    Bi_sig = [i[0] for i in Bi_sigmas]
    Bi_A = [i[0] for i in Bi_amps]
    B_ch_ave,B_ch_var = get_stats(Bi_ch)

    #-------------------------------------------------------------------------#
    # Process channel data using fit results
    #-------------------------------------------------------------------------#
    Bi_counts = fitter.get_peak_counts(Bi_ch,Bi_sig,Bi_A)
    Bi_counts = cut_outliers(Bi_counts)
    Bi_mean, Bi_var = get_stats(np.asarray(Bi_counts))

    print('K-40 <channel> = {} +/- {}'.format(K_ch_ave,K_ch_var))
    print('K-40 <N> = {} +/- {}'.format(K_mean,K_var))

    print('Bi-214 <channel> = {} +/- {}'.format(B_ch_ave,B_ch_var))
    print('Bi-214 <N> = {} +/- {}'.format(Bi_mean,Bi_var))

    #-------------------------------------------------------------------------#
    # Process weather data
    #-------------------------------------------------------------------------#
    # LBL weather station
    #location = 'KCABERKE89'
    #location = 'KCABERKE86'
    #wtimes,temps = get_weather_data(location,nhours,tstart,tstop)
    times_both,counts,temps = merge_data(bin_times,Bi_counts,wtimes,temps)
    #-------------------------------------------------------------------------#
    # Plots of everything we are interested in!
    #-------------------------------------------------------------------------#
    make_plot(bin_times,K_counts,np.sqrt(K_counts), \
              'Time','counts','K-40 counts vs Time','go','g')
    fig_name = '/Users/alihanks/Google Drive/NQUAKE_analysis/D3S/K_counts_{}_5-8.png'.format(stationID)
    plt.savefig(fig_name)

    make_plot(times,calibs,calib_err, \
              'Time','keV/channel','keV/channel vs Time','bo','b', \
              2.4,2.6)
    fig_name = '/Users/alihanks/Google Drive/NQUAKE_analysis/D3S/calibs_{}_5-8.png'.format(stationID)
    plt.savefig(fig_name)

    make_plot(bin_times,Bi_counts,np.sqrt(Bi_counts), \
             'Time','counts','Bi-214 counts vs Time','go','g')
    fig_name = '/Users/alihanks/Google Drive/NQUAKE_analysis/D3S/Bi_counts_{}_5-8.png'.format(stationID)
    plt.savefig(fig_name)

    #make_plot(Ktimes,K_ch,K_ch_errs, \
    #           'Time','1460 center channel','1460 channel vs Time','ro','r')

    #make_plot(times,Bi_ch,Bi_ch_errs, \
    #         'Time','609 center channel','609 channel vs Time','ro','r', \
    #         B_ch_ave-10*B_ch_var,B_ch_ave+10*B_ch_var)

    make_plot(temps,counts,np.sqrt(counts), \
              'Temp (F)','Bi-214 counts','Bi-214 counts vs Temp (F)','ro','r')
    fig_name = '/Users/alihanks/Google Drive/NQUAKE_analysis/D3S/Bi_counts_vs_T_{}_5-8.png'.format(stationID)
    plt.savefig(fig_name)

    plt.show()

if __name__ == '__main__':
    url = 'https://radwatch.berkeley.edu/sites/default/files/dosenet/lbl_outside_d3s.csv'
    #url = 'https://radwatch.berkeley.edu/sites/default/files/dosenet/etch_roof_d3s.csv'
    start = '2017-6-6'
    stop = '2017-5-31'
    rows = import_csv(url,start,stop)
    # number of days to look at and hours to integrate for each data point
    nhours = 1
    main(rows,nhours,start,stop)
