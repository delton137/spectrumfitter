#-----------------------------------------------------------------------
# EpsLT Peak finder - opens index of refraction data, finds peaks in
# longitudinal and transverse dielectric functions and places a line there
# 2015 D.C. Elton
#-----------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

#-----------------------------------------------------------------------
# Function to find the maxima of a dataset by looking where the sign of the slope changes
#-----------------------------------------------------------------------
def find_peaks(dataset,omegas):
    #Smooth dataset
    smoothing_length = 5
    dataset = np.convolve(dataset, np.ones(smoothing_length)/smoothing_length)
    
   # lowess = sm.nonparametric.lowess(dataset, omegas, frac=0.5)
    #dataset = lowess[:,1]

    npoints = len(dataset)
    data_shift = np.r_[0, dataset]
    diff = data_shift[0:npoints] - dataset
    peaks = []
    npoints = len(omegas) -1 
    for i in range(0, npoints):
        if np.sign(diff[i]) != np.sign(diff[i+1]):
            #check if in useful range
            if ( 580 <= omegas[i] <= 1000 ) | ( 3000 <= omegas[i] <= 3500 ): #(1500 <= omegas[i] <= 1700)
                peaks =  peaks + [(omegas[i] + omegas[i+1])/2]
    return peaks

#-----------------------------------------------------------------------
# --------------------  Main ------------------------------------------
#-----------------------------------------------------------------------

cspeed = 3*10**10
#.RI = refractive index data .DI = Dielectric data "Ice 190 K Clapp.RI"
filenames = ["Ice 200 K.RI","Ice 235 K.RI","Water 240 K.RI", "Water 253 K.RI", "Water 263 K.RI", "Water 273 K.RI","Water 300 K.DI"]
num_temps = len(filenames)
split1 = []
split2 = []
split3 = []
peak_posT = [0 for i in range(num_temps)]
peak_posL = [0 for i in range(num_temps)]
plt.close('all')

plt.subplots(num_temps,sharex=True, sharey=True)
f,  axes = plt.subplots(num_temps,sharex=True, sharey=True)

for t in range(0,num_temps):
    title = filenames[t][0:-3]
       
    #Read in data
    filename = filenames[t]
   
    #find Long. and Trans. dielectric functions
    if filename[len(filename)-2:len(filename)] == "RI": 
        data = np.loadtxt(fname=filename)
        omegas = data[:,0]
        n = data[:,1]
        k = data[:,2]
        rp = n**2 + k**2
        cp = 2*n*k
    else:
        data = np.loadtxt(fname=filename, delimiter=',')
        omegas = data[:,0]/cspeed
        rp = data[:,1]
        cp = -data[:,2]
    
    T = cp
    L = cp/(rp**2 + cp**2)
   
    #Get peak positions
    if (t == 0): 
        peak_posT[t] = [800, 3227]
        peak_posL[t] = [919, 3413]
    elif (0 < t < 6):    
        peak_posT[t] = find_peaks(T,omegas)
        peak_posL[t] = find_peaks(L,omegas)
    elif (t == 6):
        peak_posT[t] = [546,3375]
        peak_posL[t] = [750,3438]


    print(filenames[t] + " T " + str(peak_posT[t]))
    print(filenames[t] + " L " + str(peak_posL[t]))

    miny = 0 
    maxy = 2.8
    
    
    #---------------- Make the plot-----------------------
    axes[t].plot(omegas,T,'r',label=r"$\chi_T''(\omega)$")
    axes[t].plot(omegas,L,'g',label=r"$\chi_L''(\omega)$")
    axes[t].text(25, 1.75, title,fontsize=12)
    
    #add lines to plot at peak positions
    for pp in peak_posT[t]:
        axes[t].plot([pp,pp], [miny, maxy], 'r-', lw=2)
        
    for pp in peak_posL[t]:
        axes[t].plot([pp, pp], [miny, maxy], 'g-', lw=2)
   
    #Store information about LO-TO splitting
    if t <= 5:
        split1 = split1 + [peak_posL[t][0] - peak_posT[t][0]]
        split2 = split2 + [peak_posL[t][1] - peak_posT[t][1]]
       # split3 = split3 + [peak_posL[t][2] - peak_posT[t][2]]

#----------------------------------------------------------------
#--------------------------finish figure ----------------------
#----------------------------------------------------------------
#legend 
axes[3].legend(fontsize='14',frameon=False, loc='center',ncol=2,columnspacing=None)

#make subplots close to each other
f.subplots_adjust(hspace=0)

#hide x ticks for all but bottom plot.
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.setp([a.get_xticklabels() for a in f.axes], fontsize=15)
plt.setp([a.get_yticklabels() for a in f.axes], fontsize=14)

#final touches
plt.xlim([0,4000])
plt.ylim([0,2.48])
plt.yticks(np.arange(0,2, 1)) #tick spacing
plt.xlabel(r"$\omega$ (cm$^{-1}$)",fontsize=17)

plt.show() 

#----------------------------------------------------------------
#--------------------------Analyze LO-TO splittings ------------
#----------------------------------------------------------------

print("Librational LOTO Splittings: " + str(split1))
print("Bending LOTO Splittings: " + str(split2))
print("Stretching LOTO Splittings: " + str(split3))

#
for t in range(num_temps-2):
    if len(peak_posT[t]) != len(peak_posL[t]): 
        print("ERROR - missing a peak")
    else:
        LSTratio = [peak_posL[t][i]**2/peak_posT[t][i]**2 for i in [1,2,3]]
        print(filenames[t] + " " + str(LSTratio))
        print("generalized LST:" + " " +  str(LSTratio[0]*LSTratio[1]*LSTratio[2]))
            