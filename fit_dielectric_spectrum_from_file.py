from pylab import *
from scipy import optimize 
from numpy import *
from spectrum_fitter import *


def plot_model(model,dataX,dataY,handle):
    """displays a pretty plot of the model and data using matplotlib
    
    args: 
        model: a spectral_model object
        dataX: a numpy array giving the experimental x-data
        dataY: a numpy array giving the experimental y-data
        handle: an integer giving the plot window number 
    """
    figure(handle)
    plotmin = .1
    plotmax = max_freq_to_analyze
    plotomegas = linspace(10*plotmin, plotmax, 10000)

    loglog(dataX, dataY, "ro", plotomegas, model.eval(plotomegas),'g')

    #plot all of the components
    for lineshape in model.lineshapes: 
        loglog(plotomegas, lineshape(plotomegas),'g--')
    #ax.annotate('local max', xy=(3, 1),  xycoords='data')

    #xlim([plotmin,plotmax])
    #ylim([.001,max(dataY)])
#----------------------------------------------------------------------------------
#---------------------- Load dielectric data -------------------------------------
#---------------------------------------------------------------------------------- 

data = loadtxt(fname='Siegelstein.RI')
#data = loadtxt(fname='eps_omega_flex_F_TTM3F.dat')

#cspeed = 3*10**10
#data[:,0] = data[:,0]/cspeed #conv data to cm^-1
max_freq_to_analyze = 4500 #4500
maxw   = len(data[data[:,0] < max_freq_to_analyze,0])
omegas = data[0:maxw,0]

n      = data[0:maxw,1]
k      = data[0:maxw,2]
rp = n**2 + k**2
cp = 2*n*k

#rp      = data[0:maxw,1]
#cp      = -data[0:maxw,2]
  
dataX = omegas 
dataY = cp/(rp**2 + cp**2)

#----------------------------------------------------------------------------------
#---------------------- Define fit function & parameters -------------------------
#----------------------------------------------------------------------------------
#longitudinal model
print "Fitting longitudinal model..."
modelL = spectralmodel()
modelL.add(Debye([.2 , 10 ]  ,[(.001,2) , (2 ,20)]  ,  "Debye"))
#modelL.add(Debye([.01, 100]  ,[(.001,2) , (20 ,300)],"2nd Debye"))
modelL.add(DHO([2,222,200]   ,[(0,5)   ,(100 ,350)  ,(0 ,400) ],"H-bond str"))
modelL.add(DHO([.5,450,100]  ,[(.01,2) ,(380 ,650)  ,(.1,400) ],"L1"))
modelL.add(DHO([.3,760,200]  ,[(.01,2) ,(650,820)  ,(.1,1000)],"L2"))
modelL.add(DHO([1,1500,100]  ,[(.01,2) ,(1300,1700) ,(.1,500) ],"v2"))
modelL.add(DHO([1.4,3500,100],[(.01,2) ,(3000,4000) ,(.1,500) ],"v1+v3"))

modelL.fit_model(dataX,dataY)
modelL.print_model()

plot_model(modelL,dataX,dataY,1) 

##transverse model
dataY = cp
print "Fitting transverse model..."
modelT = spectralmodel()
modelT.add(Debye([71,  .5]   ,[(69,73)   ,(.1,.8)  ],"Debye"))
#modelT.add(Debye([1,   10]   ,[(.01,2)   ,(1,15)   ],"2nd Debye"))
modelT.add(DHO([1.7,200,20]  ,[(1,5)  ,(100,300)  ,(0,400)  ],"H-bond str"))
modelT.add(DHO([.5,500,300]  ,[(.01,2),(380,600)  ,(.1,400) ],"L1"))
modelT.add(DHO([.3,660,244]  ,[(.01,2),(600,800)  ,(.1,1000)],"L2"))
modelT.add(DHO([1,1500,100]  ,[(.01,2),(1300,1700),(.1,500) ],"v2"))
modelT.add(DHO([1.4,3500,100],[(.01,2),(3000,4000),(.1,500) ],"v1+v3"))
  
modelT.fit_model(dataX,dataY)
modelT.print_model()

plot_model(modelT,dataX,dataY,2) 



##----------------------------------------------------------------------------------
##---------------------- Print left hand side of LST relation ---------------------
##----------------------------------------------------------------------------------

#for lineshape in modelL.lineshapes

ratios = zeros(modelL.numlineshapes)
for i in range(modelL.numlineshapes):
	print "long.", modelL.lineshapes[i].name, "= %6.2f " % modelL.lineshapes[i].p[1]
	print "trans", modelT.lineshapes[i].name, "= %6.2f " % modelT.lineshapes[i].p[1]
	
	if modelL.lineshapes[i].type == "Debye": 
		ratios[i] = modelL.lineshapes[i].p[1]/modelT.lineshapes[i].p[1]
	
	if modelL.lineshapes[i].type == "DHO": 
		ratios[i] = (modelL.lineshapes[i].p[1]**2 + modelL.lineshapes[i].p[2]**2)/modelT.lineshapes[i].p[1]**2

set_printoptions(precision=3)

print "LST Ratios =", ratios
print "LST LHS =" , prod(ratios)
show(block=False)

