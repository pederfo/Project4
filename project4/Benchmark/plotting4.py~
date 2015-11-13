import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

"""
This is the plotting file used for output using method 5.
Note the "output.dat" on the line below indicating the datafile used.
You MUST set this yourself. For histograms there is a bins parameter where
the user sets the number of bins to divide the values between.
The value set by default is not the best for all cases but it works reasonably well.
The user is encouraged to vary this on his/her own.
"""


data = np.loadtxt("output.dat",unpack=True);
bins=abs(np.max(data)-np.min(data))
plt.figure(0)
plt.hist(data,bins,normed=1)
plt.title(r'Probability distribution of energy states')
plt.ylabel(r'P(E)',size= 14)
plt.xlabel(r'Energy',size=14)
plt.legend(loc='best')
plt.xlim([np.min(data),np.max(data)])
plt.show()
