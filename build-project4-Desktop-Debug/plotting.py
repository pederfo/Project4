from matplotlib.pyplot import *
from numpy import *
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

#Part 1

data = loadtxt("output.dat",unpack=True);

figure(0)
plot(data[0],data[2])
title(r'$<E>$ vs Monte Carlo cycles')
ylabel(r'$<E>$',size= 14)
xlabel(r'Monte Carlo cycles',size=14)
#ylim([-3,3])
#show()

figure(1)
plot(data[0],data[-1])
title(r'$<|M|>$ vs Monte Carlo cycles')
ylabel(r'$<|M|>$',size= 14)
xlabel(r'Monte Carlo cycles',size=14)
#ylim([-3,3])
show()
