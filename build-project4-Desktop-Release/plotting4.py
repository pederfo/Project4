import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

"""This file plots the percentage of accepted moves as a function of temperature
Note that the file must be named acceptedmoves.dat and you have to update the variables 
mcs, and spin with the appropriate numbers. 
Column 0 is temperature and column 1 is the number of accepted moves.
"""

data = np.loadtxt("acceptedmoves.dat",unpack=True);
mcs = 1000000
spins = 20
totalmoves = spins*spins*mcs
plt.figure(0)
plt.plot(data[0],data[1]/(totalmoves)*100)
plt.title(r'Probability distribution of energy states')
plt.ylabel(r'Percent accepted moves',size= 14)
plt.xlabel(r'Temperature [kT/J]',size=14)
plt.legend(loc='best')
plt.xlim([np.min(data[0]),np.max(data[0])])
plt.show()
