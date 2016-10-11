'''Example of how to use spectrum_fitter.py, to fit the dielectric spectrum of water.'''
__author__  = "Daniel C. Elton"
__maintainer__ = "Daniel C. Elton"
__copyright__ = "Copyright 2015, Daniel C. Elton"
__license__ = "MIT"
__status__ = "Development"

from scipy import optimize 
from numpy import *
from spectrum_fitter import *
import pickle 

#----------------------------------------------------------------------------------
#---------------------- Load dielectric data -------------------------------------
#---------------------------------------------------------------------------------- 

data = loadtxt(fname='Siegelstein.RI')
#data  = loadtxt(fname='Water 240 K.RI')

#data = loadtxt(fname='eps_omega_flex_F_TTM3F.dat')
#data = loadtxt(fname='../eps_omega_TTM3F_350K.dat')
#data = loadtxt(fname='../eps_omega_TTM3F_300K.dat')
#data = loadtxt(fname='../TTM3F_final_data/eps_omega_flex_1_300.dat')

#cspeed = 3*10**10
#data[:,0] = data[:,0]/cspeed #conv data to cm^-1

max_freq  = 1200
maxw      = len(data[data[:,0] < max_freq,0])
rawomegas = data[0:maxw,0]
min_freq  = min(rawomegas)
mid_freq  = 18
omegas    = concatenate((logspace(log10(min_freq),log10(mid_freq),30),linspace(mid_freq,max_freq,80)))

n     = data[0:maxw,1]
k     = data[0:maxw,2]
rawrp = n**2 + k**2
rawcp = 2*n*k

rp = interp(omegas,rawomegas,rawrp)
cp = interp(omegas,rawomegas,rawcp)

##----------------------------------------------------------------------------------
##---------------------- Define fit function & parameters -------------------------
##----------------------------------------------------------------------------------
##transverse model
modelT = spectralmodel()
modelT.add(Debye([71,  .5]   ,[(1,80)   ,(.4,.8)  ],"Debye"))
modelT.add(Debye([71,  .5]   ,[(1,80)   ,(.4,.8)  ],"Debye"))
modelT.add(Debye([2,   6.44]   ,[(.01,10)   ,(1,15)   ],"2nd Debye"))
#modelT.add(StretchedExp([71, 8, 1]   ,[(0,20000),(5,15), (0,1) ],"StretchedExp"))
#modelT.add(StretchedExp([200, 6.44, .91],[(0,200000),(1,15), (0,1) ],"2nd StretchedExp"))
#modelT.add(PowerLawDebye([71, .5, 1, 1.2],[(0,100),(.3,.7), (0,100), (1,2)],"PowLawDebye"))
modelT.add(BrendelDHO([1, 50, 10, 2],[(0,100  ),(1,75),(1,300),(.1,100)],"Brendel"))
modelT.add(BrendelDHO([1, 165, 10,50],[(0,100  ),(75,200),(1,300),(1,100)],"Brendel"))
modelT.add(BrendelDHO([.3,460,100,50],[(.01,100),(400,520),(1,1000),(1,300)],"Brendel L2"))
modelT.add(BrendelDHO([.3,650,100,50],[(.01,100),(520,690),(1,1000),(1,300)],"Brendel L2"))
#modelT.add(BrendelDHO([.3,680,244,50],[(.01,100),(650,720),(1,1000),(1,300)],"Brendel L3"))
#modelT.add(Debye([2,   1]   ,[(.01,4)   ,(.5,15)   ],"2nd Debye"))
#modelT.add(Debye([2,   30]   ,[(.01,4)   ,(1,100)   ],"3rd Debye"))
#modelT.add(DHO([2,60 ,200]   ,[(0,5)   ,(10  ,100)  ,(1 ,400) ],"H-bond bend"))
#modelT.add(DHO([2,150,150]   ,[(.1,4)   ,(150 ,250)  ,(10 ,900) ],"H-bond str."))
#modelT.add(DHO([.5,500,300]  ,[(.01,2),(380,550)  ,(.1,400) ],"L1"))
#modelT.add(DHO([.3,600,100]  ,[(.01,2),(550,800)  ,(1,1000)],"L2"))
#modelT.add(DHO([.3,680,244]  ,[(.01,2),(650,720)  ,(1,1000)],"L3"))
#modelT.add(DHO([1,1600,100]  ,[(.0001,.1),(1500,1700),(.1,500) ],"v2"))
#modelT.add(DHO([.21,2120,100],[(.0001,.01) ,(2000,2200) ,(10,600) ],"L+v2"))
#modelT.add(DHO([1.4,3500,100],[(.01,2),(3000,4000),(.1,500) ],"v1+v3"))
#modelT.add(constant([2]       ,[(1,11)],"eps inf"))


print("Fitting transverse model...")
modelT.fit_model(omegas,rp,cp)
 
#Optional pickling of models (save models)
#pickle.dump(modelL, open('modelL.pkl', 'wb'))
#pickle.dump(modelT, open('modelT.pkl', 'wb'))


#modelL = pickle.load(open('modelL.pkl', 'rb'))
#modelT = pickle.load(open('modelT.pkl', 'rb'))
#modelT = pickle.load(open('modelT1Debye3Brendel.pkl', 'rb'))
#modelT = pickle.load(open('modelT2Debye3Brendel.pkl', 'rb'))
#modelT = pickle.load(open('modelT1Debye4Brendel.pkl', 'rb'))
#modelT = pickle.load(open('modelT1PowLawDebye3DHO.pkl', 'rb'))

plot_model(modelT,omegas,rp,cp,4,.001,max_freq,xscale='log') 


modelT.print_model()
set_printoptions(precision=2)
print("")
print("Eps(0) = ", rp[0])


plt.show(block=True)




