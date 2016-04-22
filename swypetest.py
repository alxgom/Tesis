# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 16:20:03 2016

@author: Alexis
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tik#
from scipy.integrate import odeint
#from comparefase import comparembfase
from intmbfase import intmbfase as integ
#from intfaserev import intfaserev as intrev
#from intmbfase2 import intmbfase2 as integ2
from scipy.signal import argrelextrema
from time import localtime

#%%
#Read parameters from 'Params.txt', set.
txtfile=open("Params.txt")
tempvar=txtfile.readlines()

'''parameters for normalization'''
a= float(tempvar[1].split()[1])
gperp= float(tempvar[2].split()[1])  #gamma perpendicular, loss rate
scale=1*(10.**6)/gperp #scale to micro seconds
wscale=1000*gperp/(10.**6)#scale frequency to khz

'''parameters for the equation'''
kr= float(tempvar[3].split()[1])
k=kr/gperp #normalized loss rate
mu= float(tempvar[4].split()[1])
Dphi0= float(tempvar[5].split()[1]) #phase shift [-pi,pi]
d=1.0 #detuning
gr= float(tempvar[6].split()[1])
g=gr/gperp #*((2*pi)**2) #sigma parallel, normalized loss rate
D0=a*k/mu #Poblation
m= float(tempvar[7].split()[1]) #modulation amplitud [0,1]
wf= float(tempvar[8].split()[1])

'''parameters to compare with the results'''
w_res=np.sqrt(k*g*((D0*mu/k)-1.))*wscale #resonance frequency
a=D0*mu/k
w=np.sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale #Relaxation oscilations frequency
wf_real=wf*wscale

txtfile.close()#deberia cerrarlo aca?

'''Swype parameters'''
wfmax=0.00420
wfmin=0.00380
h=0.00001

txtfile=open("Status.txt",'r')
tempvar=txtfile.readlines()
stat= int(tempvar[0].split()[1])
print 'Status=', stat
txtfile.close()#deberia cerrarlo aca?

#%%
init=0
def  initial(init, time, y):
    intime=500.*15*10**(-6)*gperp #integration time FOR TRANSITORY
    if init==0:
        '''User defined initial condition'''
        timeinit = np.arange(0., intime, 50.)
        dfxinit=[1., 1.] 
        dfyinit=[1.,  -1.9]  
        drxinit=[1.,   1.]
        dryinit=[1.,  -1.9] 
        ddeltainit=[6.65973518e+03]
        yinit=np.array(dfxinit+dfyinit+drxinit+dryinit+ddeltainit)
    if init==1:
        '''initial condition from last simulation'''
        timeinit = np.arange(time[-1] ,intime*5/15+time[-1] , 15.)
        yinit=y[-1]
    return yinit, timeinit

def editline(file,n_line,text):
    with open(file) as infile:
        lines = infile.readlines()
    lines[n_line] = text+' \n'
    with open(file, 'w') as outfile:
        outfile.writelines(lines)        

#%%        
def swipe(m,k,mu,Dphi0,d,g,D0,wf):
    if stat==0:    
        strobo_map=[[0], [0]]
        time=np.array([0, 0])
        y=np.array([[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0]])
        yinit, time=initial(0, time, y)#first run
        y, time=integ(yinit,time,k,mu,Dphi0,d,g,D0,m,wf)
        wf=wfmax
    elif stat==1:  
        y=np.fromfile('yinitials.in',dtype=np.float64)
        time=np.fromfile('tinitials.in',dtype=np.float64)
        strobo_map=np.fromfile('strobo_map.in',dtype=np.float64)  
        
    with open('Status.txt','w') as f:
        f.write('Stat: 1')
   
    '''swipe''' 
    while wf>wfmin:
        #error print
        if wf<wf-h: print 'man, wf esta aumentando, se hace infinito el loop.' 
        wf_real=wf*wscale
        print wf
        '''initial conditions'''
        yinit, time=initial(1, time, y)
        y, time=integ(yinit,time,k,mu,Dphi0,d,g,D0,m,wf)
        yinit, time=initial(1, time, y)
        '''integration'''
        y, time=integ(yinit,time,k,mu,Dphi0,d,g,D0,m,wf)
        '''intensitys'''
        intensity_ex=np.sqrt(y[:,0]**2+y[:,1]**2)
        intensity_ey=np.sqrt(y[:,2]**2+y[:,3]**2)
        intensity=np.sqrt(y[:,0]**2+y[:,1]**2+y[:,2]**2+y[:,3]**2)
        '''map'''
        #map_time=(time)+(2*pi/wf)
        #peak_coor=argrelextrema(np.cos(wf*time), np.greater)#find peaks index
        peak_max=list(set(intensity[argrelextrema(np.cos(wf*time), np.greater)[0]]))#intensity strobo. puedo scar el set
        w_peaks=list(wf*wscale*np.ones_like(peak_max))#vector or m, the same lenght as peak_max
        strobo_map[0]=strobo_map[0]+w_peaks
        strobo_map[1]=strobo_map[1]+peak_max
        wf=wf-h   #ste next wf 

        editline('Params.txt',8,'wf: %f' %wf)#update next wf to file
        
        binwrite=open('strobo_map.in','wb')
        np.ndarray.tofile(np.array(strobo_map),'strobo_map.in')
        binwrite.close()

        binwrite=open('yinitals.in','wb')
        np.ndarray.tofile(y[-1],'yinitials.in')
        binwrite.close()
        
        binwrite=open('tinitals.in','wb')
        np.ndarray.tofile(np.array(time[-1]),'tinitials.in')
        binwrite.close()        
    return strobo_map  
    
    

#%%    
strobo_map=swipe(m,k,mu,Dphi0,d,g,D0,wf)    
with open('Status.txt','w') as f:
        f.write('Status: 0') #ends the run, sets status back to 0.


#%%    
'''Plots!''' 
save=False #set True if i want to save files automatically
pi=np.pi

#%%    
fig5=plt.figure()
fig5.suptitle('Phase spase Map. Max intensity vs. Wf', fontsize=12, fontweight='bold')
ax1 = fig5.add_subplot(111)
ax1.plot(np.array(strobo_map[0])*wscale,strobo_map[1],'.b')
ax1.set_xlabel('w [kHz]')
ax1.set_ylabel('Strobo Intensity |E|')
ax1.set_xlim(np.min(strobo_map[0])*wscale, np.max(strobo_map[0])*wscale)
ax1.set_ylim(0, 2.5)
plt.text(-0.1,-.32, "\n Parameters: $m= $ %s , $\Delta \phi_0=$ %s ,  $\Omega=$ %.2f khz \n $k$=%.2f khz, $\mu'=$ %s , $\delta= $ %.2e , $\gamma_{||}=$ %s , $D_0=$ %s, $A=$ %.1f " % (m, Dphi0, w_res ,k,mu, d, g, D0, a), fontsize=11, transform=ax1.transAxes)   
plt.subplots_adjust(bottom=0.22)
fig5.set_size_inches(7,7)
if save==True: 
    fname='%d_%d_%d-%d.%d.%d-max_vs_w.png' % localtime()[0:6]
    fig5.savefig(fname, dpi = 100)# when saving, specify the DPI
plt.show
