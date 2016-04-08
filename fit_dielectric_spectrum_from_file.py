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

max_freq_to_analyze = 1000
maxw   = len(data[data[:,0] < max_freq_to_analyze,0])
omegas = data[0:maxw,0]

n      = data[0:maxw,1]
k      = data[0:maxw,2]
rp = n**2 + k**2
cp = 2*n*k

Ldatarp = rp/(rp**2 + cp**2)
Ldatacp = cp/(rp**2 + cp**2)

##----------------------------------------------------------------------------------
##---------------------- Define fit function & parameters -------------------------
##----------------------------------------------------------------------------------
##longitudinal model
modelL = spectralmodel()
#modelL.add(Debye([.2 ,  .4 ]  ,[(.0005,1) , (.1 ,20)]  ,  "low freq Debye"))
modelL.add(Debye([.001 , 10 ]   ,[(.001,1)  , (.5 ,50)]   ,"Debye"))
modelL.add(Debye([5 , 100 ]   ,[(.001,1)  , (.5 ,50)]   ,  "2nd Debye"))
modelL.add(Debye([.1, 100]   ,[(.0001,1) , (20 ,200)],     "3rd Debye"))
modelL.add(DHO([2,60 ,60 ]   ,[(0,5)   ,(10  ,100)  ,(1  ,400) ],"H-bond bend"))
modelL.add(DHO([2,222,200]   ,[(0,5)   ,(200 ,300)  ,(10 ,900) ],"H-bond str."))
modelL.add(DHO([.5,450,100]  ,[(.01,2) ,(380 ,700)  ,(.1,400) ],"L1"))
modelL.add(DHO([.3,660,244]  ,[(.01,2) ,(300,700)  ,(1,1000)],"L2"))
#modelL.add(DHO([.3,760,200]  ,[(.01,2) ,(650,800)  ,(1,1000)],"L3"))
modelL.add(DHO([1,1500,100]  ,[(.0001,.1) ,(1300,1700) ,(10,600) ],"v2"))
modelL.add(DHO([.21,2120,100],[(.0001,.01) ,(2000,2200) ,(10,600) ],"L+v2"))
modelL.add(DHO([1.4,3500,100],[(.01,2) ,(3000,4000) ,(.1,500) ],"v1+v3"))

##transverse model
modelT = spectralmodel()
modelT.add(Debye([71,  .5]   ,[(68,73)   ,(.3,.6)  ],"Debye"))
modelT.add(Debye([2,   10]   ,[(.01,4)   ,(1,15)   ],"2nd Debye"))
modelT.add(Debye([2,   30]   ,[(.01,4)   ,(1,100)   ],"3rd Debye"))
modelT.add(DHO([2,60 ,200]   ,[(0,5)   ,(10  ,100)  ,(1 ,400) ],"H-bond bend"))
modelT.add(DHO([2,150,150]   ,[(.1,4)   ,(150 ,250)  ,(10 ,900) ],"H-bond str."))
modelT.add(DHO([.5,500,300]  ,[(.01,2),(380,700)  ,(.1,400) ],"L1"))
modelT.add(DHO([.3,600,100]  ,[(.01,2),(500,650)  ,(1,1000)],"L2"))
#modelT.add(DHO([.3,680,244]  ,[(.01,2),(650,720)  ,(1,1000)],"L3"))
modelT.add(DHO([1,1600,100]  ,[(.0001,.1),(1500,1700),(.1,500) ],"v2"))
modelT.add(DHO([.21,2120,100],[(.0001,.01) ,(2000,2200) ,(10,600) ],"L+v2"))
modelT.add(DHO([1.4,3500,100],[(.01,2),(3000,4000),(.1,500) ],"v1+v3"))
  
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

plot_model(modelL,omegas,Ldatacp,1,.001,max_freq_to_analyze) 
plot_model(modelT,omegas,cp,2,.001,max_freq_to_analyze)
plot_model(modelL,omegas,Ldatacp,3,.001,max_freq_to_analyze,scale='log') 
plot_model(modelT,omegas,cp,4,.001,max_freq_to_analyze,scale='log') 


##----------------------------------------------------------------------------------
##----- Printout all frequencies in system and left side of gLST equation  --------
##----------------------------------------------------------------------------------
ratios = zeros(modelL.numlineshapes)
sumL = 0
sumT = 0 
print("      name      |   f    |       freq       |      tau (ps) ")

if (modelL.numlineshapes == modelT.numlineshapes): 
	for i in range(modelL.numlineshapes):
		if modelL.lineshapes[i].type == "Debye":
			print(("long. %11s %6.2f  %6.2f (%5.2f ps)" % (modelL.lineshapes[i].name, modelL.lineshapes[i].p[0], modelL.lineshapes[i].p[1], 33.34/(2*3.141*modelL.lineshapes[i].p[1])  )))
		
		if modelT.lineshapes[i].type == "Debye":
			print(("trans %11s %6.2f  %6.2f (%5.2f ps)" % (modelT.lineshapes[i].name, modelT.lineshapes[i].p[0], modelT.lineshapes[i].p[1], 33.34/(2*3.141*modelT.lineshapes[i].p[1]) )))
			ratios[i] = modelL.lineshapes[i].p[1]/modelT.lineshapes[i].p[1]		
	
		if (modelL.lineshapes[i].type == "DHO") or (modelL.lineshapes[i].type == "VanVleck"): 
			print(("long. %11s %6.2f  %6.2f + %6.2f i  %6.3f  " % (modelL.lineshapes[i].name, modelL.lineshapes[i].p[0], modelL.lineshapes[i].p[1], modelL.lineshapes[i].p[2],  33.34/(2*3.141*modelL.lineshapes[i].p[2]))))
	
		if (modelT.lineshapes[i].type == "DHO") or (modelT.lineshapes[i].type == "VanVleck"): 
			print(("trans %11s %6.2f  %6.2f + %6.2f i  %6.3f " % (modelT.lineshapes[i].name, modelT.lineshapes[i].p[0], modelT.lineshapes[i].p[1], modelT.lineshapes[i].p[2],  33.34/(2*3.141*modelT.lineshapes[i].p[2]))))	
			ratios[i] = (modelL.lineshapes[i].p[1]**2 + modelL.lineshapes[i].p[2]**2)/modelT.lineshapes[i].p[1]**2

	
		sumL = sumL + modelL.lineshapes[i].p[0]
		sumT = sumT + modelT.lineshapes[i].p[0]
else: 
	print("ERROR: number of lineshapes in Longitudinal model not equal to number in Transverse model")
		
print("")
set_printoptions(precision=2)
print(("LST Ratios =", ratios))
print(("LST LHS = %6.2f" % prod(ratios)))
print("")
print(("Sum of tran. oscillator strengths = %6.2f" % sumT)) 
print(("Sum of long. oscillator strengths = %6.2f" % sumL))
print("")
print(("long. RMS error = %6.3f" % modelL.RMS_error))
print(("trans RMS error = %6.3f" % modelT.RMS_error))

show(block=True)


