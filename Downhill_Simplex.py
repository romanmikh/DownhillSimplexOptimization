# see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin.html
# Function to minimize (adapted to a power law)
def func(pars, x, ydat, dy):
    exponent,const=pars
    ymod = const*x**exponent
    return sum((ydat-ymod)**2/dy**2)

# Function to minimize (adapted to a power law for converging ratio)
def funcratio(pars, x, ydat, dy):
    exponent,const=pars
    ymod = const*x**exponent + 1  # +1 to ensure convergence for large values of the primary axis
    return sum((ydat-ymod)**2/dy**2)


# Now fit a line to the data with fmin by minimizing chi^2
from scipy.optimize import fmin
#x = y = err = []
#fit = fmin(func,[1,1],args=(x, y, err))
#yfit = fit[1]*x**fit[0]
#print(" Input: exponent = %5.2f; const = %5.2f"%(1, 1))
#print(" Fit:   exponent = %5.2f; const = %5.2f"%(fit[0], fit[1]))



"""-----------------------------EXAMPLE-----------------------------"""

# code used to produce part of fig. 7.1
# make sure 'r2_Runs.csv' is in the appropriate subdirectory

from pylab import sqrt
import pylab as plt
import numpy as np
import csv

# Import our and Sambaran's merger data
f = open('r2_Runs.csv')
csv_f = csv.reader(f)

banmerg = []
ourmerg = []
for row in csv_f:
    banmerg.append(row[6])
    ourmerg.append(row[11])

banmerg = np.array(list(np.float_(banmerg[2:]))) # 2 bc skipping column title and the first run with 0 mergers
ourmerg = np.array(list(np.float_(ourmerg[2:])))


# Input Sambaran's Data and our data
x = np.array((15000, 20000, 30000, 50000, 70000, 75000, 100000)) # masses

# Ratio and error propagation
errban = np.sqrt(banmerg)
ct = np.array((1,5,26,23,12,2,4)) # number of points in each mass bin for Sambaran's data
ctcum = np.array((1,6,32,55,67,69,73))

# sum of Banerjee's data for bins of masses
summ = np.zeros(len(ct))
for i in range(len(ct)-1):
    summ[i+1] = np.sum(banmerg[ctcum[i]:ctcum[i+1]])
summ[0] = banmerg[0]

avgban = np.zeros(len(ct))
for i in range(len(ct)):
    avgban[i] = summ[i]/ct[i] # averages, first part of error propagation equation


# sum of errors for bins of masses
sumerr = np.zeros(len(ct))
for i in range(len(ct)-1):
    sumerr[i+1] = np.sum(errban[ctcum[i]:ctcum[i+1]])
sumerr[0] = errban[0]

banerr = np.zeros(len(ct))
for i in range(len(ct)):
    banerr[i] = sumerr[i]/ct[i] # average errors



# Our data points, errors are stdevs, since no intrinsic errors exist
oursum = np.zeros(len(ct))
for i in range(len(ct)-1):
    oursum[i+1] = np.sum(ourmerg[ctcum[i]:ctcum[i+1]])
oursum[0] = ourmerg[0]

avgour = np.zeros(len(ct))
for i in range(len(ct)):
    avgour[i] = oursum[i]/ct[i]

# our errors
ourerr = np.zeros(len(ct))
for i in range(len(ct)-1):
    ourerr[i+1] = np.std(ourmerg[ctcum[i]:ctcum[i+1]])
ourerr[0] = np.std(ourmerg[0])



# Now fit a line to the data with fmin by minimizing chi^2
fit = fmin(func,[1,1],args=(x,avgban,banerr))
yfit = fit[1]*x**fit[0]
print(" Input: exponent = %5.2f; const = %5.2f"%(1, 0))
print(" Fit:   exponent = %5.2f; const = %5.2f"%(fit[0], fit[1]))
