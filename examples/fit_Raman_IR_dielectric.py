'''Example of how to use spectrum_fitter.py, to fit three spectra (Raman, transverse dielectric function, longitudinal dielectric function) and make a plot of all three'''
__author__  = "Daniel C. Elton"
__maintainer__ = "Daniel C. Elton"
__copyright__ = "Copyright 2015, Daniel C. Elton"
__license__ = "MIT"
__status__ = "Development"

from pylab import *
from scipy import optimize 
from numpy import *
from spectrumfitter.spectrumfitter import *
from spectrumfitter.spectralmodel import *
import pickle 

#---------------------- Load data -------------------------------------
eps_data = loadtxt(fname='water_full_Siegelstein.RI')

Raman_data = loadtxt(fname='water_300K_Raman_Castner.dat')


#----------------process dielectric data ------------------------------
max_freq  = 1200
maxw      = len(eps_data[eps_data[:,0] < max_freq,0])
rawomegas = eps_data[0:maxw,0]
min_freq  = min(rawomegas)
mid_freq  = 18
omegas    = rawomegas #concatenate((logspace(log10(min_freq),log10(mid_freq),30),linspace(mid_freq,max_freq,80)))
n     = eps_data[0:maxw,1]
k     = eps_data[0:maxw,2]
rawrp = n**2 + k**2
rawcp = 2*n*k

rp = rawrp #interp(omegas,rawomegas,rawrp)
cp = rawcp #interp(omegas,rawomegas,rawcp)

denom = rp**2 + cp**2 	
(Lrp,Lcp) = (rp/denom, cp/denom)

#----------------process Raman data ------------------------------------
maxw = len(Raman_data[Raman_data[:,0] < max_freq,0])
Raman_omegas = Raman_data[0:maxw,0]
Raman = Raman_data[0:maxw,1]/5000


##------------ Define fit functions & parameters -----------------------
##transverse model
modelT = SpectralModel()
modelT.add(Debye([71,  .5]   ,[(1,80)   ,(.4,.8)  ],"Debye"))
#modelT.add(Debye([71,  .5]   ,[(1,80)   ,(.4,.8)  ],"Debye"))
modelT.add(Debye([2,   6.44]   ,[(.01,10)   ,(1,15)   ],"2nd Debye"))
#modelT.add(StretchedExp([71, 8, 1]   ,[(0,20000),(5,15), (0,1) ],"StretchedExp"))
#modelT.add(StretchedExp([200, 6.44, .91],[(0,200000),(1,15), (0,1) ],"2nd StretchedExp"))
#modelT.add(PowerLawDebye([71, .5, 1, 1.2],[(0,100),(.3,.7), (0,100), (1,2)],"PowLawDebye"))
#modelT.add(BrendelDHO([1, 50, 10, 2],[(0,100  ),(1,75),(1,300),(.1,100)],"Brendel"))
modelT.add(BrendelDHO([1, 165, 10,50],[(0,100  ),(75,200),(1,300),(1,100)],"Brendel"))
#modelT.add(BrendelDHO([.3,460,100,50],[(.01,100),(400,520),(1,1000),(1,300)],"Brendel L2"))
#modelT.add(BrendelDHO([.3,650,100,50],[(.01,100),(520,750),(1,1000),(1,300)],"Brendel L2"))
modelT.add(DHO([.3,460,100],[(.01,100),(400,600),(1,1000)],"DHO L1"))
modelT.add(DHO([.3,650,100],[(.01,100),(520,750),(1,1000)],"DHO L2"))
#modelT.add(BrendelDHO([.3,680,244,50],[(.01,100),(650,720),(1,1000),(1,300)],"Brendel L3"))
#modelT.add(Debye([2,   1]   ,[(.01,4)   ,(.5,15)   ],"2nd Debye"))
#modelT.add(Debye([2,   30]   ,[(.01,4)   ,(1,100)   ],"3rd Debye"))
#modelT.add(DHO([2,60 ,200]   ,[(0,5)   ,(10  ,100)  ,(1 ,400) ],"H-bond bend"))
#modelT.add(DHO([2,150,150]   ,[(.1,4)   ,(150 ,250)  ,(10 ,900) ],"H-bond str."))
modelT.add(constant([2]       ,[(1,11)],"eps inf"))

##longitudinal model
modelL = SpectralModel()
modelL.add(Debye([1 , 10 ]   ,[(.001,1)  , (.5 ,50)]   ,"Debye"))
modelL.add(Debye([.1 , 10 ]   ,[(0,2)  , (.5 ,50)]   ,  "2nd Debye"))
#modelL.add(DHO([2,60 ,60 ]   ,[(0,5)   ,(10  ,100)  ,(1  ,400) ],"H-bond bend"))
modelL.add(DHO([2  ,222,200]   ,[(0,5)   ,(100 ,300)  ,(1 ,900) ],"H-bond str."))
modelL.add(DHO([.2,450,100]  ,[(0,2) ,(380 ,600)  ,(.1,400) ],"L1"))
modelL.add(DHO([.1,660,244]  ,[(0,2) ,(600,770)  ,(1,1000)],"L2"))
#modelL.add(BrendelDHO([.3,680,244,50],[(.01,100),(650,720),(1,1000),(1,300)],"Brendel L3"))
modelL.add(constant([2]       ,[(1,10)],"eps inf"))



#Raman model
modelR = SpectralModel()
#modelR.add(DHO([.1,65,100]   ,[(0,1)   ,(10 ,100)  ,(10 ,900) ],"H-bond bend"))
#modelR.add(DHO([.1,150,100]   ,[(0,1)   ,(100 ,200)  ,(10 ,900) ],"H-bond str."))
modelR.add(DHO([.1,433,250]  ,[(0,.5),(300,500)  ,(.1,400) ],"L1"))
modelR.add(DHO([.04,660,250]  ,[(0,.5),(400,750)  ,(1,1000)],"L2"))
modelR.add(DHO([.015,802,200]  ,[(0,.5) ,(650,800)  ,(1,1000)],"L3"))
#modelR.add(constant([0]       ,[(0,.4)],"eps inf"))


print("Fitting transverse model...")
modelT.fit_model(omegas,rp,cp)

print("Fitting longitudinal model...")
modelL.fit_model(omegas,Lrp,Lcp)

print("Fitting Raman model...")#
#modelR.fit_model(Raman_omegas,Raman,Raman)
 
#Optional pickling of models (save models)
#pickle.dump(modelL, open('modelL.pkl', 'wb'))
#pickle.dump(modelT, open('modelT.pkl', 'wb'))
#pickle.dump(modelR, open('modelR.pkl', 'wb'))


#modelL = pickle.load(open('modelL.pkl', 'rb'))
#modelT = pickle.load(open('modelT2Debye3Brendel.pkl', 'rb'))
#modelT = pickle.load(open('2Debye2DHO3Libr3DHOT.pkl', 'rb'))

#plot_model(modelT,omegas,rp,cp,1,xmin=min_freq,xmax=max_freq,xscale='log',yscale='log')
#plot_model(modelT,omegas,rp,cp,2,xmin=2,xmax=max_freq,ymin=-.1,ymax=6,xscale='linear',yscale='linear') 
#plot_model(modelL,omegas,Lrp,Lcp,1,xmin=min_freq,xmax=max_freq,xscale='log',yscale='log')
#plot_model(modelL,omegas,Lrp,Lcp,2,xmin=2,xmax=max_freq,ymin=-.1,ymax=6,xscale='linear',yscale='linear') 

#plot_model(modelR,Raman_omegas,Raman,Raman,1,xmin=min_freq,xmax=max_freq,xscale='log',yscale='log')
#plot_model(modelR,Raman_omegas,Raman,Raman,2,xmin=2,xmax=max_freq,xscale='linear',yscale='linear') 


f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

(rpfit, cpfit) = modelT(omegas)
ax1.plot(omegas,cp,'r',omegas,cpfit)
for lineshape in modelT.lineshapes: 
	(rpPart, cpPart) = lineshape(omegas)
	ax1.plot(omegas, cpPart ,'g')

ax1.set_ylim([0,3])

ax2.plot(omegas,Lcp,'r',linewidth=3)
for lineshape in modelL.lineshapes: 
	(Lrp, Lcp) = lineshape(omegas)
	ax2.plot(omegas, Lcp ,'g')

ax3.plot(Raman_omegas,Raman,'r',linewidth=3)
for lineshape in modelR.lineshapes: 
	(Rrp, Rcp) = lineshape(omegas)
	ax3.plot(omegas, Rcp ,'g')

ax3.set_ylim([0,.4])
ax3.set_ylim([0,.4])
ax3.tick_params(axis='x',labelsize=20)

plt.xlabel(r"$\omega$ (cm$^{-1}$)",fontsize=25)

#ax1.text(500, 1.35, 'Transverse dielectric spectra',fontsize=17)
#ax2.text(500, 1.35, 'Longitudinal dielectric spectra',fontsize=17) 
#ax3.text(500, 1.35, 'Raman spectra',fontsize=17)
#plt.xlim([10,4000])
#plt.xlim([-.5,2])

modelL.print_model()
modelT.print_model()
#modelR.print_model()
#

set_printoptions(precision=2)
print("")
print("Eps(0) = ", rp[0])


show(block=True)




