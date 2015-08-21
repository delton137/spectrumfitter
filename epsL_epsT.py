#-----------------------------------------------------------------------
# Make plot comparing three sources of epsL, epsT data
# 2015 D.C. Elton
#-----------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

cspeed = 3*10**10
miny=0
maxy=2.4

data = np.loadtxt(fname='eps_omega_flex_F_TIP4P2005f.dat')
omegas1 = data[:,0]/cspeed
rp = data[:,1]
cp = -data[:,2]
TTip4p = cp
LTip4p = cp/(rp**2 + cp**2)
   
data = np.loadtxt(fname='eps_omega_TTM3F.dat')
omegas2 = data[:,0]/cspeed
rp = data[:,1]
cp = -data[:,2]
TTTM3F = cp
LTTM3F = cp/(rp**2 + cp**2)
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

data = np.loadtxt(fname='eps_omega_combined.csv',delimiter=",")
omegas3 = data[:,0]/cspeed
rp = data[:,1]
cp = -data[:,2]
Texpt = rp
Lexpt = cp/(rp**2 + cp**2)

#-------------------- Three subplots sharing both x/y axes -----------------
plt.close('all')
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
#plt.ylim([miny, maxy])
ax1.plot(omegas1,TTip4p,'r',label=r"$\chi_T''(\omega)$")
ax1.plot(omegas1,LTip4p,'g',label=r"$\chi_L''(\omega)$")
ax1.legend(fontsize='x-large',frameon=False, loc='upper right')
ax2.plot(omegas2,TTTM3F,'r',omegas2,LTTM3F,'g')
ax3.plot(omegas3,Texpt,'r',omegas3,Lexpt,'g')
ax1.text(1500, 1.35, 'TIP4P/2005f',fontsize=17)
ax2.text(1500, 1.35, 'TTM3F',fontsize=17) 
ax3.text(1500, 1.35, 'Experimental',fontsize=17)

plt.xlim([10,4000])
plt.xlim([-.5,2])


#Add lines 
peak_posT = [0, 0, 0]
peak_posL = [0, 0, 0]
peak_posT[0] = [626, 3320]
peak_posL[0] = [879, 3320]
peak_posT[1] = [583, 3430]
peak_posL[1] = [782, 3645]
peak_posT[2] = [550, 3363]
peak_posL[2] = [750, 3449]


for pp in peak_posT[0]:
  ax1.plot([pp,pp], [miny, maxy], 'r-', lw=2)
for pp in peak_posT[1]:
  ax2.plot([pp,pp], [miny, maxy], 'r-', lw=2)
for pp in peak_posT[2]:
  ax3.plot([pp,pp], [miny, maxy], 'r-', lw=2)
for pp in peak_posL[0]:
  ax1.plot([pp,pp], [miny, maxy], 'g-', lw=2)
for pp in peak_posL[1]:
  ax2.plot([pp,pp], [miny, maxy], 'g-', lw=2)
for pp in peak_posL[2]:
  ax3.plot([pp,pp], [miny, maxy], 'g-', lw=2)



# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.setp([a.get_xticklabels() for a in f.axes], fontsize=15)
plt.setp([a.get_yticklabels() for a in f.axes], fontsize=15)
plt.xlabel(r"$\omega$ (cm$^{-1}$)",fontsize=17)



'''
plt.subplot(3, 1, 1)
plt.plot(omegas1,TTip4p,'r',omegas1,LTip4p,'g')
plt.xlabel(r'$n(\omega)$',fontsize=16)
plt.ylabel(r'Im$\lbrace\varepsilon(\omega)\rbrace$', fontsize=16)
plt.xlim([0,4000])
plt.ylim([0,2])
plt.title('TIP4P/2005f')

plt.subplot(3, 1, 2)
plt.plot(omegas2,TTTM3F,'r',omegas2,LTTM3F,'g')
plt.ylabel(r'Im$\lbrace\varepsilon(\omega)\rbrace$', fontsize=16)
plt.xlim([0,4000])
plt.ylim([0,2])
plt.title('TTM3F')

plt.subplot(3, 1, 3)
plt.plot(omegas3,Texpt,'r',omegas3,Lexpt,'g')
plt.ylabel(r'Im$\lbrace\varepsilon(\omega)\rbrace$', fontsize=16)
plt.xlim([0,4000])
plt.ylim([0,2])
plt.title('expt')
'''

plt.show()