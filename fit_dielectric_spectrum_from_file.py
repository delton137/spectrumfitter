from pylab import *
from scipy import optimize 
from numpy import *
from spectrum_fitter import *
import pickle 

#----------------------------------------------------------------------------------
#---------------------- Load dielectric data -------------------------------------
#---------------------------------------------------------------------------------- 

data = loadtxt(fname='Siegelstein.RI')
#data = loadtxt(fname='eps_omega_flex_F_TTM3F.dat')

#cspeed = 3*10**10
#data[:,0] = data[:,0]/cspeed #conv data to cm^-1
max_freq_to_analyze = 4250 #4500
maxw   = len(data[data[:,0] < max_freq_to_analyze,0])
omegas = data[0:maxw,0]

n      = data[0:maxw,1]
k      = data[0:maxw,2]
rp = n**2 + k**2
cp = 2*n*k
show(block=False)

#rp      = data[0:maxw,1]
#cp      = -data[0:maxw,2]
Ldata = cp/(rp**2 + cp**2)
Tdata = cp

##----------------------------------------------------------------------------------
##---------------------- Define fit function & parameters -------------------------
##----------------------------------------------------------------------------------
#longitudinal model
modelL = spectralmodel()
#modelL.add(BRO([.2 , 10, 10 ]  ,[(.001,2) , (.01,20), (.1 ,20)]  ,  "Debye"))
modelL.add(Debye([.2 , 10 ]  ,[(.001,2) , (2 ,20)]  ,  "Debye"))
modelL.add(Debye([.1, 100]   ,[(.001,2) , (20 ,200)],"2nd Debye"))
modelL.add(DHO([2,60 ,60 ]   ,[(0,5)   ,(10  ,100)  ,(1  ,400) ],"H-bond bend"))
modelL.add(DHO([2,222,200]   ,[(0,5)   ,(100 ,300)  ,(10 ,500) ],"H-bond str."))
modelL.add(DHO([.5,450,100]  ,[(.01,2) ,(380 ,650)  ,(.1,400) ],"L1"))
modelL.add(DHO([.3,760,200]  ,[(.01,2) ,(650,820)  ,(.1,1000)],"L2"))
modelL.add(DHO([1,1500,100]  ,[(.01,2) ,(1300,1700) ,(.1,500) ],"v2"))
modelL.add(DHO([1.4,3500,100],[(.01,2) ,(3000,4000) ,(.1,500) ],"v1+v3"))

##transverse model
modelT = spectralmodel()
modelT.add(Debye([71,  .5]   ,[(69,73)   ,(.1,.8)  ],"Debye"))
modelT.add(Debye([1,   10]   ,[(.01,2)   ,(1,15)   ],"2nd Debye"))
modelT.add(DHO([2,60 ,200]   ,[(0,5)   ,(10  ,100)  ,(1 ,400) ],"H-bond bend"))
modelT.add(DHO([2,222,200]   ,[(0,5)   ,(100 ,300)  ,(10 ,500) ],"H-bond str."))
modelT.add(DHO([.5,500,300]  ,[(.01,2),(380,600)  ,(.1,400) ],"L1"))
modelT.add(DHO([.3,660,244]  ,[(.01,2),(500,700)  ,(.1,1000)],"L2"))
modelT.add(DHO([1,1500,100]  ,[(.01,2),(1300,1700),(.1,500) ],"v2"))
modelT.add(DHO([1.4,3500,100],[(.01,2),(3000,4000),(.1,500) ],"v1+v3"))
  
print "Fitting longitudinal model..."
modelL.fit_model(omegas,Ldata)

print "Fitting transverse model..."
modelT.fit_model(omegas,Tdata)
 
#1 Debye + 5DHO
modelL.setparams([0.14716527009116862, 16.663114929754016, 0.13977221495572345, 266.26785634863182, 324.96617318843602, 0.083976259671020057, 634.79034989125341, 349.24686688652497, 0.055620926702680422, 791.00660256746517, 203.16930952582692, 0.01, 1649.4477783763575, 184.83666777080384, 0.018893044599961276, 3449.0976363271166, 259.79569131498704])
modelT.setparams([73.0, 0.64998712207262099, 2.1554588510176402, 190.94392910415954, 396.01620705563181, 0.38740036835447894, 539.75636626673759, 281.18778513438099, 0.17987654276372583, 690.67335496511498, 244.01062804336505, 0.013223803766142527, 1644.9290701596349, 78.405182960132578, 0.061509563470691704, 3368.4991788840525, 265.93348589445338])
modelL.print_model()
modelT.print_model()

#Optional pickling of models (save models)
pickle.dump(modelL, open('modelL.pkl', 'wb'))
pickle.dump(modelT, open('modelT.pkl', 'wb'))


#modelL = pickle.load(open('modelL.pkl', 'rb'))
#modelT = pickle.load(open('modelT.pkl', 'rb'))

plot_model(modelL,omegas,Ldata,1,.001,max_freq_to_analyze) 
plot_model(modelT,omegas,Tdata,2,.001,max_freq_to_analyze) 



##----------------------------------------------------------------------------------
##---------------------- Print left hand side of LST relation ---------------------
##----------------------------------------------------------------------------------

#for lineshape in modelL.lineshapes

ratios = zeros(modelL.numlineshapes)
sumL = 0
sumT = 0 
print "      name      |   f    |       freq"
for i in range(modelL.numlineshapes):
	if modelL.lineshapes[i].type == "Debye":
		print "long. %11s %6.2f  %6.2f" % (modelL.lineshapes[i].name, modelL.lineshapes[i].p[0], modelL.lineshapes[i].p[1])
		print "trans %11s %6.2f  %6.2f" % (modelT.lineshapes[i].name, modelT.lineshapes[i].p[0], modelT.lineshapes[i].p[1])
		ratios[i] = modelL.lineshapes[i].p[1]/modelT.lineshapes[i].p[1]
	
	if (modelL.lineshapes[i].type == "DHO") or (modelL.lineshapes[i].type == "BRO"): 
		print "long. %11s %6.2f  %6.2f + %6.2f i " % (modelL.lineshapes[i].name, modelL.lineshapes[i].p[0], modelL.lineshapes[i].p[1], modelL.lineshapes[i].p[2])
		print "trans %11s %6.2f  %6.2f + %6.2f i " % (modelT.lineshapes[i].name, modelT.lineshapes[i].p[0], modelT.lineshapes[i].p[1], modelT.lineshapes[i].p[2])

		ratios[i] = (modelL.lineshapes[i].p[1]**2 + modelL.lineshapes[i].p[2]**2)/modelT.lineshapes[i].p[1]**2
		
	sumL = sumL + modelL.lineshapes[i].p[0]
	sumT = sumT + modelT.lineshapes[i].p[0]
		
print ""
set_printoptions(precision=2)
print "LST Ratios =", ratios
print "LST LHS = %6.2f" % prod(ratios)
print ""
print "Sum of long oscillator strengths = %6.2f" % sumT 
print "Sum of tran oscillator strengths = %6.2f" % sumL
print ""
print "long. RMS error = %6.3f" % modelL.RMS_error
print "trans RMS error = %6.3f" % modelT.RMS_error

show(block=True)

