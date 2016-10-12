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

#maxw   = len(data[data[:,0] < max_freq_to_analyze,0])
#omegas = data[0:maxw,0]

max_freq_to_analyze = 1300
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
##longitudinal model
modelL = SpectralModel()
modelL.add(Debye([.1,  10]  ,[(.01,.5) , (.1 ,50)] , "Debye"))
modelL.add(Debye([.01 , 20 ] ,[(.0001,.3) , (.1 ,150)] , "2nd Debye"))
modelL.add(Debye([.01 , 150] ,[(.0001,.3) , (.1 ,300)], "3rd Debye"))
modelL.add(DHO([2,60 ,60]   ,[(0,5)   ,(10  ,100)  ,(1  ,400)],"Hbond Bend"))
modelL.add(DHO([.1,230,200]  ,[(.001,.4)  ,(200 ,300)  ,(10 ,900) ],"Hbond Str"))
modelL.add(DHO([.1,450,100]  ,[(.001,.4) ,(380 ,700)  ,(.1,400) ]," L1"))
modelL.add(DHO([.1,750,244]  ,[(.001,.4) ,(450,800)  ,(1,1000)]," L2"))
modelL.add(DHO([.3,760,200]  ,[(.01,2) ,(650,800)  ,(1,1000)]," L3"))
#modelL.add(BrendelDHO([2,60 ,60,10 ]   ,[(0,5)   ,(10  ,100)  ,(1  ,400), (1,100)],"Brendel Hbond Bend"))
#modelL.add(BrendelDHO([.1,230,200,10]  ,[(.001,.4)  ,(200 ,300)  ,(10 ,900),(1,100) ],"Brendel Hbond Str"))
#modelL.add(BrendelDHO([.1,450,100,10]  ,[(.001,.4) ,(380 ,700)  ,(.1,400),(1,100) ],"Brendel L1"))
#modelL.add(BrendelDHO([.1,750,244,10]  ,[(.001,.4) ,(450,800)  ,(1,1000),(1,100)],"Brendel L2"))
#modelL.add(BrendelDHO([.3,760,200,10]  ,[(.01,2) ,(650,800)  ,(1,1000),(1,100)],"Brendel L3"))
#modelL.add(DHO([1,1500,100]  ,[(.0001,.1) ,(1300,1700) ,(10,600) ],"v2"))
#modelL.add(DHO([.21,2120,100],[(.0001,.01) ,(2000,2200) ,(10,600) ],"L+v2"))
##modelL.add(BrendelDHO([1.4,3500,100,20],[(.01,2) ,(3000,4000) ,(.1,500),(1,20) ],"v1+v3"))
#modelL.add(DHO([1.4,3500,100],[(.01,2) ,(3000,4000) ,(.1,500) ],"v1+v3"))
modelL.add(constant([1]       ,[(.01,1.1)],"eps inf L"))


##transverse model
modelT = SpectralModel()
modelT.add(Debye([69,  .55]   ,[(65,73)   ,(.3,.65)    ],"Debye"))
modelT.add(Debye([2,   2]    ,[(.001,10)   ,(.1,10)   ],"2nd Debye"))
modelT.add(Debye([2,   5]   ,[(.001,5)   ,(3,50)   ],"3rd Debye"))
modelT.add(DHO([2,60 ,200]   ,[(0,5)   ,(10  ,100)  ,(1 ,400) ],"H-bond bend"))
modelT.add(DHO([2,150,150]   ,[(.1,4)   ,(150 ,250)  ,(10 ,900) ],"H-bond str."))
modelT.add(DHO([.5,400,300]  ,[(.01,2),(340,500)  ,(.1,400) ],"L1"))
modelT.add(DHO([.3,500,100]  ,[(.01,2),(450,700)  ,(1,1000)],"L2"))
modelT.add(DHO([.3,650,100]  ,[(.01,2),(550,750)  ,(1,1000)],"L3"))
#modelT.add(BrendelDHO([1, 50, 10, 2],[(0,100  ),(1,75),(1,300),(.1,100)],"Brendel Hbond Bend"))
#modelT.add(BrendelDHO([1, 165, 10,50],[(0,100  ),(75,200),(1,300),(1,100)],"Brendel Hbond Str"))
#modelT.add(BrendelDHO([.3,460,100,40],[(.01,100),(400,520),(1,500),(1,150)],"Brendel L1"))
#modelT.add(BrendelDHO([.3,650,100,40],[(.01,100),(520,750),(1,500),(1,150)],"Brendel L2"))
#modelT.add(BrendelDHO([.3,680,244,40]  ,[(.01,2),(650,720)  ,(1,1000),(1,150)],"Brendel L3"))
#modelT.add(DHO([1,1600,100]  ,[(.0001,.1),(1500,1700),(.1,500) ],"v2"))
#modelT.add(DHO([.21,2120,100],[(.0001,.01) ,(2000,2200) ,(10,600) ],"L+v2"))
##modelT.add(BrendelDHO([1.4,3500,200,20],[(.01,2),(3000,4500),(.1,500),(1,100) ],"Brendelv1+v3"))
#modelT.add(DHO([1.4,3500,200],[(.01,2),(3000,4500),(.1,500)],"v1+v3"))
modelT.add(constant([2]       ,[(1,11)],"eps inf"))


##----------------------------------------------------------------------------------
##---------------------- Fitting the models----------------------------------------
##----------------------------------------------------------------------------------

# ----------- fitting longitudinal and transverse models seperately ---------------
print("Fitting longitudinal model...")
modelL.fit_model(omegas,Ldatarp,Ldatacp)

print("Fitting transverse model...")
modelT.fit_model(omegas,rp,cp)

 
#Optional pickling of models (save models)
#pickle.dump(modelL, open('modelL.pkl', 'wb'))
#pickle.dump(modelT, open('modelT.pkl', 'wb'))

#pickle.dump(modelL, open('TTM3FL.pkl', 'wb'))
#pickle.dump(modelT, open('TTM3FL.pkl', 'wb'))

#modelL = pickle.load(open('modelL.pkl', 'rb'))
#modelT = pickle.load(open('modelT.pkl', 'rb'))

#modelL = pickle.load(open('TTM3FL.pkl', 'rb'))
#modelT = pickle.load(open('TTM3FL.pkl', 'rb'))

#modelL = pickle.load(open('2DebyeHstr2Lib3DHOL.pkl', 'rb'))
#modelT = pickle.load(open('2DebyeHstr2Lib3DHOT.pkl', 'rb'))

#-------------------- fit both models together with gLST constraint -------------------

#print("doing gLST constrained fit")
#fit_model_gLST_constraint(modelL, modelT, omegas, rp, cp)

#-------------------- priting out parameters and plotting ------------------------------
print_gLST_LHS_stuff(modelL, modelT)

plot_model(modelL,omegas,Ldatarp,Ldatacp,1,.001,max_freq_to_analyze,title='Longitudinal eps(omega)') 
plot_model(modelT,omegas,rp,cp,2,.001,max_freq_to_analyze,title='Transverse eps(omega)')

plot_model(modelL,omegas,Ldatarp,Ldatacp,3,.001,max_freq_to_analyze,xscale='log',title='Longitudinal eps(omega)') 
plot_model(modelT,omegas,rp,cp,4,.001,max_freq_to_analyze,xscale='log',title='Transverse eps(omega)') 

plt.show(block=True)