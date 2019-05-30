import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

var = input("Enter the .csv filename to calibrate: ")

csv = np.genfromtxt(var, delimiter= ",").T
summed = np.sum(csv, axis=1)


def fitFunc(x, a, m, s, c):
    return a * np.exp(-(x - m)**2 / (2 * s**2)) + c
def linBgFitFunc(x, a, m, s, c, b):
    return a * np.exp(-(x - m)**2 / (2 * s**2)) + c + b * x
def find(xSlice, xshift, trymax=20, trymu=200, trysig=100, trybg=5):
    xmu = np.mean(xSlice)
    xsig = np.std(xSlice)
    xxdata = range(len(xSlice))
    trydata = fitFunc(xSlice, np.max(xSlice), xmu, xsig, np.max(xSlice) + 50)
    p0 = [trymax,trymu,trysig,trybg]
    xpopt, xpcov = curve_fit(fitFunc, xxdata, xSlice, p0)
    xchannel = xshift + int(xpopt[1])
    return xchannel
def linBgFind(xSlice, xshift, trymax=20, trymu=200, trysig=100, trybg=5, trylin=-20):
    xmu = np.mean(xSlice)
    xsig = np.std(xSlice)
    xxdata = range(len(xSlice))
    p0 = [trymax,trymu,trysig,trybg, trylin]
    xpopt, xpcov = curve_fit(linBgFitFunc, xxdata, xSlice, p0)
    xchannel = xshift + int(xpopt[1])
    return xchannel
'''
def showFindFit(xSlice, xshift, trymax=20, trymu=200, trysig=100, trybg=5, lin=1):
    xmu = np.mean(xSlice)
    xsig = np.std(xSlice)
    xxdata = range(len(xSlice))
    p0 = [trymax,trymu,trysig,trybg, lin]
    xpopt, xpcov = curve_fit(linBgFitFunc, xxdata, xSlice, p0)
    xchannel = xshift + int(xpopt[1])
    return linBgFitFunc(xxdata, *xpopt)
'''
def findChannel(shift, range):
    cSlice = summed[shift:shift+range]
    return find(cSlice, shift)
def lbFindChannel(shift, range):
    cSlice = summed[shift:shift+range]
    return linBgFind(cSlice, shift, 1200, 60, 80, 20)
Bi2 = findChannel(1600, 300)
Bi4 = findChannel(900, 200)
Bi5 = findChannel(540, 100)
Pb = lbFindChannel(250, 130)
Th = lbFindChannel(60, 150)

channels = [Bi2, Bi4, Bi5, Pb, Th]
energies = [1120.3, 609.3, 352, 185.7, 92.6]

def polyfit(x, m, b, r):
    return r * x*x + m*x + b
p0 = [1, .6, 2]
xpopt, xpcov = curve_fit(polyfit, channels, energies, p0)
print("xpopt:", xpopt)
plt.plot(polyfit(range(max(channels)), *xpopt))
plt.plot(channels, energies, 'r.')
plt.show()
print("best fit line:", xpopt[2], "x^2 + ", xpopt[0], "x + ", xpopt[1])
print("Channels of peaks: ", channels)