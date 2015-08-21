from pylab import *
from scipy import optimize

 
 
data = loadtxt(fname='Siegelstein.RI')
data[:,0] = 10000/data[:,0]#convert um into cm^-1
data = flipud(data)
savetxt('Siegelstein.RI2',data,'%10.6f,%10.6f,%10.6f' )
