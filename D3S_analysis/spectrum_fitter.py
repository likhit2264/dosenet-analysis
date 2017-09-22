import csv
import argparse
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

# allow for reloading libraries in Python 2/3
try:
    reload
except NameError:
    from importlib import reload

import spectra_fitting_tools as fitter
reload(fitter)

#--------------------------------------------------------------------------#
# Process input data
#--------------------------------------------------------------------------#
def make_int(lst): 
    '''
    Makes all entries of a list an integer
    '''
    y = []
    for i in lst:
        i = i.translate({ord(c): None for c in '[]'})
        y.append(int(i))
    return y

def make_array(lst,lower,upper=5000): 
    '''
    Makes list into an array. Also splices out the irrelevant stuff 
    for a spectra
    '''
    y = np.asarray(make_int(lst[lower:upper]))
    return y

def import_local_csv(file,lower,upper=5000):
    data = []
    with open(file) as f:
        reader = csv.reader(f)
        rows = [r for r in reader]

    for row in rows:
        data.append(make_array(row,lower,upper))
    return data

def main(input_file,limits,plot=False):
    data = import_local_csv(input_file,limits[0],limits[1])
    integrated = sum(data)
    if plot:
        fig = plt.figure()
        fig.patch.set_facecolor('white')
        plt.title('Spectra')
        plt.xlabel('channels')
        plt.ylabel('counts')
        plt.xlim(limits[0],limits[1])
        x = list(range(limits[0],len(integrated)+limits[0]))
        plt.plot(x,integrated,'b:',label='data')
        plt.yscale('log')
        plt.show()

    else:
        mean,sigma,amp = fitter.single_peak_fit(integrated,limits[0],limits[1],20.0,50,True)
        print('mean = {}, sigma = {}, amp = {}'.format(mean,sigma,amp))
        peak_counts = fitter.get_peak_counts(mean[0],sigma[0],amp[0])
        raw_counts = fitter.get_gross_counts(integrated,0,5000)
        print('peak counts = {}'.format(peak_counts))
        print('gross counts = {}'.format(raw_counts))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--input-file", type=str, dest='input_file', default='data.csv')
    parser.add_argument('--range', nargs='+', type=int, dest='range', default=[0,5000])
    parser.add_argument('--plot', dest='plot', action='store_true', default=False)
    info = parser.parse_args()

    main(info.input_file,info.range,info.plot)
