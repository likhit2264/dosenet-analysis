import csv
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

#import weather_data_tools as weather
#reload(weather)
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
        y.append(int(i))
    return y

def make_array(lst,lower,upper=5000): 
    '''
    Makes list into an array. Also splices out the irrelevant stuff 
    for a spectra
    '''
    y = np.asarray(make_int(lst[lower:upper]))
    return y

def import_local_csv(file):
    with open(file) as f:
        reader = csv.reader(f)
        rows = [r for r in reader]

    data = make_array(rows,12)
    return data

def main(input_file,range):
    data = import_local_csv(input_file)
    integrated = sum(data)
    mean,sigma,amp = fitter.single_peak_fit(integrated,range[0],range[1],5.0)
    peak_counts = fitter.get_peak_counts(mean[0],sigma[0],amp[0])
    raw_counts = get_gross_counts(integrated,range[0],range[1])

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--input-file", type = string, dest='input_file', default='data.csv')
    parser.add_argument('--range', nargs='+', type=int, dest='range', default=[0,5000])
    info = parser.parse_args()

    main(info.input_file,info.range)
