from pylab import *
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

max_freq_to_analyze = 100
maxw   = len(data[data[:,0] < max_freq_to_analyze,0])
omegas = data[0:maxw,0]

n      = data[0:maxw,1]
k      = data[0:maxw,2]
rp = n**2 + k**2
cp = 2*n*k

#rp      = data[0:maxw,1]
#cp      = -data[0:maxw,2]

Ldata = cp/(rp**2 + cp**2)
Tdata = cp

##----------------------------------------------------------------------------------
##---------------------- Define fit function & parameters -------------------------
##----------------------------------------------------------------------------------
##transverse model
modelT = spectralmodel()
modelT.add(Debye([71,  .5]   ,[(68,73)   ,(.3,.7)  ],"Debye"))
#modelT.add(StretchedExp([71, 8, 1]   ,[(0,20000),(5,15), (0,1) ],"StretchedExp"))
#modelT.add(PowerLawDebye([71, .5, 1, 1],[(0,100),(.3,.7), (0,100), (0,2)],"Power law Debye"))
modelT.add(BrendelDHO([4, 1, 1, .5],[(0,80),(.1,100),(.01,100),(.5,.5)],"Brendel"))
#modelT.add(Debye([2,   1]   ,[(.01,4)   ,(.5,15)   ],"2nd Debye"))
#modelT.add(Debye([2,   30]   ,[(.01,4)   ,(1,100)   ],"3rd Debye"))
#modelT.add(DHO([2,60 ,200]   ,[(0,5)   ,(10  ,100)  ,(1 ,400) ],"H-bond bend"))
#modelT.add(DHO([2,150,150]   ,[(.1,4)   ,(150 ,250)  ,(10 ,900) ],"H-bond str."))
#modelT.add(DHO([.5,500,300]  ,[(.01,2),(380,700)  ,(.1,400) ],"L1"))
#modelT.add(DHO([.3,600,100]  ,[(.01,2),(500,650)  ,(1,1000)],"L2"))
#modelT.add(DHO([.3,680,244]  ,[(.01,2),(650,720)  ,(1,1000)],"L3"))
#modelT.add(DHO([1,1600,100]  ,[(.0001,.1),(1500,1700),(.1,500) ],"v2"))
#modelT.add(DHO([.21,2120,100],[(.0001,.01) ,(2000,2200) ,(10,600) ],"L+v2"))
#modelT.add(DHO([1.4,3500,100],[(.01,2),(3000,4000),(.1,500) ],"v1+v3"))
modelT.add(constant([2]       ,[(1,11)],"eps inf"))


print "Fitting transverse model..."
modelT.fit_model(omegas,rp,cp)
 
#Optional pickling of models (save models)
#pickle.dump(modelL, open('modelL.pkl', 'wb'))
#pickle.dump(modelT, open('modelT.pkl', 'wb'))


#modelL = pickle.load(open('modelL.pkl', 'rb'))
#modelT = pickle.load(open('modelT.pkl', 'rb'))

#plot_model(modelL,omegas,Ldata,3,.001,max_freq_to_analyze,scale='log') 
plot_model(modelT,omegas,rp,cp,4,.001,max_freq_to_analyze,scale='log') 


#Ldata = cp/(rp**2 + cp**2)

##----------------------------------------------------------------------------------
##----- Printout all frequencies in system and left side of gLST equation  --------
##----------------------------------------------------------------------------------
sumT = 0 
print "      name      |   f    |       freq       |      tau (ps) "

for i in range(modelT.numlineshapes):
	if modelT.lineshapes[i].type == "Debye":
		print "trans %11s %6.2f  %6.2f (%5.2f ps)" % (modelT.lineshapes[i].name, modelT.lineshapes[i].p[0], modelT.lineshapes[i].p[1], 33.34/(2*3.141*modelT.lineshapes[i].p[1]) )

	elif (modelT.lineshapes[i].type == "DHO") or (modelT.lineshapes[i].type == "BRO") or (modelT.lineshapes[i].type == "StretchedExp"): 
		print "trans %11s %6.2f  %6.2f + %6.2f i  %6.3f " % (modelT.lineshapes[i].name, modelT.lineshapes[i].p[0], modelT.lineshapes[i].p[1], modelT.lineshapes[i].p[2],  33.34/(2*3.141*modelT.lineshapes[i].p[2]))	
		
	elif modelT.lineshapes[i].type == "Constant":
		print "trans %11s %6.2f" % (modelT.lineshapes[i].name, modelT.lineshapes[i].p[0])
	
	else:
		print modelT.lineshapes[i].p

	sumT = sumT + modelT.lineshapes[i].p[0]
print ""
set_printoptions(precision=2)
print ""
print "Sum of trans. oscillator strengths = %6.2f" % sumT
print ""
print "trans RMS error = %6.3f" % modelT.RMS_error

show(block=True)


