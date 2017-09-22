import numpy as np
from scipy import optimize
from scipy import asarray as ar,exp
from scipy.integrate import quad

import matplotlib.pyplot as plt

verbose = 0
#--------------------------------------------------------------------------#
# Fit Functions
#--------------------------------------------------------------------------#
def lbound(bound,par):
    return 1e4*np.sqrt(bound-par) + 1e-3*(bound-par) if (par<bound) else 0

def ubound(bound,par):
    return 1e4*np.sqrt(par-bound) + 1e-3*(par-bound) if (par>bound) else 0

def bound(bounds,par):
    return lbound(bounds[0],par) + ubound(bounds[1],par)

def fixed(fix,par):
    return bound((fix,fix), par)

def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))+lbound(0,a)+lbound(0,sigma)+lbound(0,x0)

def expo(x,a,slope):
    return a*exp(x*slope)+lbound(0,a)+ubound(0,slope)

# p = [a1,mean,sigma,a2,shift,slope,const]
def gaus_plus_exp(x,p):
    return gaus(x,p[0],p[1],p[2])+expo(x,p[3],p[4])

# p = [a1,mean,sigma,slope,const]
def gaus_plus_line(x,p):
    return gaus(x,p[0],p[1],p[2])+p[3]*x+p[4]

def double_gaus_plus_exp(x,p):
    return gaus(x,p[0],p[1],p[2])+gaus(x,p[3],p[4],p[5])+expo(x,p[6],p[7])

def double_gaus_plus_line(x,p):
    return gaus(x,p[0],p[1],p[2])+gaus(x,p[3],p[4],p[5])+p[6]*x+p[7]

def peak_fitter(x,y,fit_function,pinit): 
    """
    Peak Finder for peak in specified range

    Args:
        x: data x values for fitting
        y: data y values for fitting
        fit_function: fit function
        pinit: inital parameters for fit function

    Returns:
        array of resulting fit parameters and array of fit errors
    """
    errfunc = lambda p, x, y: fit_function(x,p) - y 
    pfit,pcov,infodict,errmsg,success = \
        optimize.leastsq(errfunc, pinit, args=(x,y), \
            full_output=1, epsfcn=0.0001)

    if (len(y) > len(pinit)) and pcov is not None:
        s_sq = (errfunc(pfit, x, y)**2).sum()/(len(y)-len(pinit))
        pcov = pcov * s_sq
    else:
        pcov = 0

    error = [] 
    for i in range(len(pfit)):
        try:
          error.append(np.absolute(pcov[i][i])**0.5)
        except:
          error.append( 0.00 )
    pfit_leastsq = pfit
    perr_leastsq = np.array(error) 
    return pfit_leastsq, perr_leastsq 

def single_peak_fit(array,lower,upper,sigma,count_offset=1,make_plot=False,save_plot=False,plot_name=''):
    """
    Performs single gaussian + exponential background fit

    Args:
        array: full array of counts (spectra)
        lower,upper: bounds on spectra for window to fit inside
        count_offset: correction for shift from left edge of spectrum
        make_plot: flag for plotting fit result (diagnostic)

    Returns:
        list of fit parameters+errors
    """
    points = ar(range(lower,upper))
    count_list = list(array[:])
    counts = array

    nentries = len(points)
    mean = lower + (upper - lower)/2.0
    max_value = max(count_list)
    max_index = count_list.index(max_value)
    if max_index > points[0]+20:
        mean = max_index
    max_counts = array[0]
    min_counts = array[-1]
    if min_counts == 0:
        min_counts = 1
    amp = max_value - (max_counts - min_counts)/2.0
    slope = (np.log(min_counts)-np.log(max_counts))/(points[-1]-points[0])
    pinit = [amp,mean,sigma,array[0]*count_offset,slope]
    print('Initial parameters: {}'.format(pinit))
    pars,errs = peak_fitter(points,array,gaus_plus_exp,pinit)
    print('Fit parameters: {}'.format(pars))
    if make_plot:
        fig = plt.figure()
        fig.patch.set_facecolor('white')
        plt.title('Spectra integrated over a day')
        plt.xlabel('channels')
        plt.ylabel('counts')
        #plt.xlim(lower,upper)
        #plt.ylim(counts[-1]*.1,counts[0]*10)
        x = ar(range(0,len(array)))
        plt.plot(points,array,'b:',label='data')
        #pars = [ 2.95010675e+01, 1.06815654e+03, 6.94962149e+01, 3.89127957e+03, -4.64346847e-03]
        plt.plot(points,gaus_plus_exp(points,pars),'ro:',label='fit')
        plt.legend()
        plt.yscale('log')
        if save_plot:
            #'/Users/alihanks/Google Drive/NQUAKE_analysis/D3S/fit_plots/'
            fig_file = plot_name+'.pdf'
            plt.savefig(fig_file)
            plt.close()
        else:
            plt.show()

    if verbose:
        par_labels = ['norm','mean','sigma','amp','slope']
        for i in range(len(pars)):
            print('{}-{}: {} +/- {}'.format(par_labels[i],counter,pars[i],errs[i]))

    return [pars[1],errs[1]],[pars[2],errs[2]],[pars[0],errs[0]]

def double_peak_fit(array,counter,lower,upper,pindex=0,count_offset=1,make_plot=False,plot_name=''):
    """
    Performs double gaussian + exponential background fit

    Args:
        array: full array of counts (spectra)
        lower,upper: bounds on spectra for window to fit inside
        pindex: indication of which gaussian to get fit results for
        count_offset: correction for shift from left edge of spectrum
        make_plot: flag for plotting fit result (diagnostic)

    Returns:
        list of fit parameters+errors
    """
    points = ar(range(lower,upper))
    counts = ar(list(array[lower:upper]))

    nentries = len(points)
    mean = lower + (upper - lower)/2.0
    slope = (np.log(counts[-1])-np.log(counts[0]))/(points[-1]-points[0])
    pinit = [counts[0]/7.0,mean-5.0,3.0,counts[0]/7.0,mean+5.0,3.0, \
             counts[0]*count_offset,slope]
    pars,errs = peak_fitter(points,counts,double_gaus_plus_exp,pinit)
    if verbose:
        par_labels = ['norm1','mean1','sigma1','norm2','mean2','sigma2','amp','slope']
        for i in range(len(pars)):
            print('{}-{}: {} +/- {}'.format(par_labels[i],counter,pars[i],errs[i]))

    if make_plot:
        fig = plt.figure()
        fig.patch.set_facecolor('white')
        plt.title('Spectra integrated over a day')
        plt.xlabel('channels')
        plt.ylabel('counts')
        plt.xlim(lower,upper)
        plt.ylim(20,1000)
        x = ar(range(0,len(array)))
        plt.plot(x,array,'b:',label='data')
        plt.plot(x,double_gaus_plus_exp(x,pars),'ro:',label='fit')
        plt.legend()
        plt.yscale('log')
        fig_file = '/Users/alihanks/Google Drive/NQUAKE_analysis/D3S/fit_plots/'+plot_name+'_fit_'+str(counter)+'.pdf'
        plt.savefig(fig_file)
        plt.close()

    mean = [pars[1],errs[1]]
    sigma = [pars[2],errs[2]]
    amp = [pars[0],errs[0]]
    if (pindex==1 and pars[4] > pars[1]) or (pindex==0 and pars[4] < pars[1]):
        mean = [pars[4],errs[4]]
        sigma = [pars[5],errs[5]]
        amp = [pars[3],errs[3]]
        if errs[1] > errs[4]:
            mean[1] = errs[1]
        if abs(pars[2]-pars[5])/pars[2] > 0.8:
            mean[1] = 150
    return mean,sigma,amp

def get_peak_counts(mean,sigma,amp):
    count,err = quad(gaus,0,5000,args=(amp,mean,sigma))
    return count,err

def get_all_peak_counts(means,sigmas,amps):
    '''
    Calculate the area under a gaussian curve (estimate of counts in that peak)

    Arguments:
      - list of guassian means
      - list of guassian widths
      - list of gaussian amplitudes

    Returns:
      - list of counts from resulting gaussian integrations 
    '''
    counts = []
    for i in range(len(means)):
        count,err = get_peak_counts(means[i],sigmas[i],amps[i])
        counts.append(count)
    return counts

def get_gross_counts(array,lower,upper):
    counts = sum(array[lower:upper])
    return counts

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
    # print(rows)

    datatz = rows[-1][1].tzinfo
    date_itr = tstart
    times = []
    means = []
    sigmas = []
    amps = []
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
            mean,sigma,amp = fit_function(integrated,*fit_args)
            means.append(mean)
            sigmas.append(sigma)
            amps.append(amp)
            times.append(integration[int(len(integration)/2)][1])

    return times,means,sigmas,amps
