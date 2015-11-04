from matplotlib.pyplot import *
from numpy import *
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

#Part 1

data = loadtxt("output.dat",unpack=True);

figure(0)
plot(data[0],data[2],label=r'$<E>$')
plot(data[0],data[-1],label=r'$<|M|>$')
title(r'Expectation values vs Monte Carlo cycles')
ylabel(r'Expectation value',size= 14)
xlabel(r'Monte Carlo cycles',size=14)
legend(loc='best')
ylim([-2.5,1.5])
#show()

C_V = data[3]/data[1]*data[1]
chi = data[5]/data[1]
figure(1)
plot(data[0],C_V,label=r'$C_V$')
plot(data[0],chi,label=r'$\chi$')
legend(loc='best')
title(r'Magnetic susceptibility and Specific heat',size=14)
ylabel(r'$C_v$ and $\chi$',size=14)
xlabel(r'Monte Carlo cycles',size=14)
show()
