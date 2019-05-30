import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
testnum = int(input("Enter test number (1, 2): "))
detnum = int(input("Enter det number (1, 8): "))
gentext = "RPi_data/Test_" + str(testnum) + "_p1_g" + str(detnum) + "_2019-05-28_D3S.csv"

csv = np.genfromtxt(gentext, delimiter= ",").T
summed = np.sum(csv, axis=1)

'''
This is a helper fuction that looks at each index and checks if it is a peak.
REMOVED:::::Does not look at values under 1/4 of np.average(data):::::::
'''
def checkShape(i, data, r, e):
    sweep = [data[i + dx] for dx in range(-r, r+1)]
    prev=sweep[r]
    if not prev == max(sweep):# or prev < np.average(data)/4:
        return False
#     if not prev > np.average(sweep) * 1.5:
#         return False
    e = e * 2
    # ^because the code checks r indices to the left and right
    for k in range(1, r+1):
        if e < 0:
            #print(e)
            return False
        if sweep[r-k] > prev:
            e = e - 1
        prev = sweep[r-k]
    prev=sweep[r]
    for k in range(1, r+1):
        if e < 0:
            return False
        if sweep[r+k] > prev:
            e = e - 1
        prev = sweep[r+k]
    return e >= 0

'''
Takes in a summed peak count, a peak range, and an error allowance and returns possible peaks.
Peak range is the number of values the function will look at on either side
Error allowance is the number of values within the peak range that are allowed to not fit a downwards slope
'''
def sweepLeft(data, r, e):
    peaks = []
    index = r
    while index < len(data) - r:
        if checkShape(index, data, r, e):
            peaks.append(index)
            index = index + r - e//2
        else:
            index += 1
    return peaks
peakRange = int(input("Enter a peak range: "))
errAllo = int(input("Enter an error allowance: "))


ldots = sweepLeft(summed, peakRange, errAllo)
print("returned peaks:", ldots)
print("len peaklist:", len(ldots))
#print(len(ldots))
#print(np.average(summed)/4)
x=np.arange(len(summed))
plt.plot(summed)
#plt.plot(x, np.average(summed)/4 + 0*x)
plt.plot(ldots, summed[ldots], 'ro')
plt.yscale('log')
plt.show()