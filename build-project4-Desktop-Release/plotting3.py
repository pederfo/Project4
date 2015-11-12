import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

"""This file plots a histogram for the probability
distribution of the energy states
This file needs an output file with one single column 
with the energy for a given monte carlo cycle plotted on each line.
THIS ONLY WORKS FOR ONE SINGLE TEMPERATURE AT A TIME
"""
data = np.loadtxt("output.dat",unpack=True);
bins=abs(np.max(data)-np.min(data))
plt.figure(0)
#plt.hist(data,bins,normed=1)
plt.hist(data,bins,normed=1)
plt.title(r'Probability distribution of energy states')
plt.ylabel(r'P(E)',size= 14)
plt.xlabel(r'Energy',size=14)
plt.legend(loc='best')
plt.xlim([np.min(data),np.max(data)])
plt.show()
