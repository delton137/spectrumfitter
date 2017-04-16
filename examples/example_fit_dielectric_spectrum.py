__author__  = "Daniel C. Elton"
__maintainer__ = "Daniel C. Elton"
__copyright__ = "Copyright 2015, Daniel C. Elton"
__license__ = "MIT"
__status__ = "Development"

from scipy import optimize 
from numpy import *
from spectrumfitter.spectrumfitter import *
from spectrumfitter.spectralmodel import *
import pickle 
#----------------------------------------------------------------------------------
#---------------------- Load dielectric data -------------------------------------
#---------------------------------------------------------------------------------- 


data = loadtxt(fname='water_full_Siegelstein.RI')

max_freq_to_analyze = 2000

maxw      = len(data[data[:,0] < max_freq_to_analyze,0])
rawomegas = data[0:maxw,0]
min_freq  = min(rawomegas)
mid_freq  = 10
omegas    = concatenate((logspace(log10(min_freq),log10(mid_freq),150),linspace(mid_freq,max_freq_to_analyze,140)))

n      = data[0:maxw,1]
k      = data[0:maxw,2]

rawrp = n**2 + k**2
rawcp = 2*n*k

rp = interp(omegas,rawomegas,rawrp)
cp = interp(omegas,rawomegas,rawcp)

Ldatarp = 1.0 - rp/(rp**2 + cp**2)
Ldatacp = cp/(rp**2 + cp**2)

##----------------------------------------------------------------------------------
##---------------------- Define fit function & parameters -------------------------
##----------------------------------------------------------------------------------

##transverse model
modelT = SpectralModel()
modelT.add(Debye([69,  .55]   ,[(65,73)   ,(.3,.65)    ],"Debye"))
modelT.add(Debye([2,   2]    ,[(.0001,10)   ,(.5,10)   ],"2nd Debye"))
modelT.add(DHO([2,60 ,200]   ,[(0,10)   ,(10  ,100)  ,(1 ,400) ],"H-bond bend"))
#modelT.add(DHO([2,150,150]   ,[(.1,4)   ,(150 ,250)  ,(10 ,900) ],"H-bond str."))
#modelT.add(BrendelDHO([1, 50, 10, 2],[(0,100  ),(1,75),(1,300),(.1,100)],"Brendel Hbond Bend"))
#modelT.add(BrendelDHO([1, 165, 10,50],[(0,100  ),(75,200),(1,300),(1,100)],"Brendel Hbond Str"))
#modelT.add(BrendelDHO([.3,460,100,40],[(.01,100),(400,520),(1,500),(1,150)],"Brendel L1"))
#modelT.add(BrendelDHO([.3,650,100,40],[(.01,100),(520,750),(1,500),(1,150)],"Brendel L2"))
#modelT.add(DHO([1,1600,100]  ,[(.0001,.1),(1500,1700),(.1,500) ],"v2"))
#modelT.add(DHO([1,2120,100],[(.0001,.01) ,(2000,2400) ,(10,600) ],"L+v2"))
##modelT.add(BrendelDHO([1.4,3500,200,20],[(.01,2),(3000,4500),(.1,500),(1,100) ],"Brendelv1+v3"))
#modelT.add(DHO([1.4,3500,200],[(.01,2),(3000,4500),(.1,500)],"v1+v3"))
modelT.add(constant([0]       ,[(1,11)],"eps inf"))


##---------------------- Fitting the model ----------------------------------------

print("Fitting transverse model...")
 
modelT.fit_model(omegas,rp,cp)


## Optional pickling of models (save models)
#pickle.dump(modelT, open('modelT.pkl', 'wb'))

#modelT = pickle.load(open('modelT.pkl', 'rb'))


##---------------------- Plotting the model ----------------------------------------
plot_model(modelT,omegas,rp,cp,4,.05,max_freq_to_analyze,title='Transverse') 

modelT.print_model()



plt.show(block=True)