from matplotlib.pyplot import *
from numpy import *
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

"""
This is the plotting file used for output using method 2.
Note the "output.dat" on the line below indicating the datafile used.
You MUST set this yourself.
"""


data = loadtxt("L40.dat",unpack=True);
figure(0)
plot(data[0],data[1],label=r'$<E>$')
plot(data[0],data[5],label=r'$<|M|>$')
title(r'Expectation values vs Temperature')
ylabel(r'Expectation value',size= 14)
xlabel(r'Temperature[$kT/J$]',size=14)
legend(loc='best')
xlim([min(data[0]),max(data[0])])

C_V = data[2]/(data[0]*data[0])
chi = data[4]/data[0]

figure(1)
plot(data[0],C_V,label=r'$C_V$')
plot(data[0],chi,label=r'$\chi$')
title(r'Specific heat and susceptibility values vs Temperature')
ylabel(r'',size= 14)
xlabel(r'Temperature[$kT/J$]',size=14)
legend(loc='best')
xlim([min(data[0]),max(data[0])])
show()

