# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 10:23:10 2022

@author: roelvink

Python version of runup formulations Van Ormondt, 2O21 
"""

import numpy as np
from scipy import interpolate as intp

from .disper import disper, disper_fentonmckee

class runup_vo21:
    
    def __init__(self, hm0, tp, ztide, sl1, sl2, sl1opt):
                
        self.hm0= hm0
        self.tp=tp
        self.h=self.tp*10
        
        k = disper(2*np.pi/self.tp, self.h, 9.81)
        # n   = 0.5*(1+2*k*self.h/np.sinh(2*k*self.h))
        # c   = 9.81*self.tp/(2*np.pi)*np.tanh(k*self.h)
        # cg  = n*c
        
        self.l1=2*np.pi/k
        self.steepness=self.hm0/self.l1
        
        if np.size(ztide)==1:
            ztide=np.zeros(np.size(np.size(self.hm0)))+ztide
        
        if np.size(sl2)==1:
            sl2=np.zeros(np.size(self.hm0))+sl2
            
        if sl1opt == 'ad':
            gambr=1
            fh=1/gambr
            surfslope2=np.zeros(np.size(ztide))
            surfslope2[0]=(fh*self.hm0)/(fh*self.hm0/self.surfslope)**1.5
            for itide in range(np.size(ztide)):
                if ztide[itide]==0:
                    tmp=0
                else:
                    xxx= np.arange(-1000, 5005, 5)
                    yyy0=-sl1[itide]*xxx**(2/3)
                    yyy0[np.argwhere(xxx<0)]=-xxx[np.argwhere(xxx<0)]*sl2[itide]
                    yyy=yyy0-ztide[itide]
                    ii1=np.argwhere(yyy<0)[0]
                    ii2=np.argwhere(yyy<-self.hm0[itide]/gambr)[0]
                    surfslope2[itide]=(yyy[ii1]-yyy[ii2])/(xxx[ii2]-xxx[ii1])
                
        elif sl1opt == 'slope':
            surfslope2=sl1
        elif sl1opt == 'xz':
            gambr=1
            xxxx=np.arange(sl1['x'][0], sl1['x'][-1], 1)
            zzzz=intp.interp1d(sl1['x'], sl1['z'])(xxxx)
            sl1['x']=xxxx
            sl1['z']=zzzz
            surfslope2=np.zeros(np.size(ztide))
            xxx=sl1['x']
            yyy0=sl1['z']
            for itide in range(np.size(ztide)):
              
                yyy=yyy0-ztide[itide]
                ii1=np.argwhere(yyy<0)[0]
                ii2=np.argwhere(yyy<-self.hm0[itide]/gambr)[0]
                surfslope2[itide]=(yyy[ii1]-yyy[ii2])/(xxx[ii2]-xxx[ii1])
        self.sl2=sl2
        self.ztide= ztide
        self.surfslope=surfslope2           
        self.ksis=self.surfslope/(np.sqrt(self.hm0/self.l1)) 
                    
    def compute_r2p(self, drspr, phi):
       
        ksi1=self.surfslope/(np.sqrt(self.hm0/self.l1))
        ksi2=self.sl2/(np.sqrt(self.hm0/self.l1))

        self.compute_setup(drspr, phi)
        self.compute_hm0_lf(drspr, phi)
        self.compute_hm0_hf(drspr, phi)
        
        self.r2p=self.setup+0.82396*np.sqrt((0.82694*self.hm0_lf)**2+(0.73965*self.hm0_hf)**2)*ksi2**0.15201*ksi1**(-0.086635)
        
    def compute_setup(self, drspr, phi):
        beta=[None]*7
        beta[0]=4.0455506e+00 
        beta[1]=2.3740615e-02 
        beta[2]=2.0340287e+00 
        beta[3]=4.6497588e-01 
        beta[4]=7.0244541e-01 
        beta[5]=5.0000000e-01 
        beta[6]=-3.1727583e-01 
        
        ksib=self.sl2/(np.sqrt(self.hm0/self.l1))
    
        v=beta[0]*self.hm0*(beta[1]+np.exp(-beta[2]*self.ksis**beta[6]*ksib**beta[3])*ksib**beta[4]) 
        fac1=self.compute_dirspreadfac_setup(drspr) 
        fac2=self.compute_directionfac_setup(phi) 
        self.setup=v*fac1*fac2 
        
    def compute_hm0_lf(self, drspr, phi):
        beta=[None]*8
        beta[0]=3.4547125e+00
        beta[1]=5.8790748e-01
        beta[2]=3.6906975e+00
        beta[3]=2.3378556e-01
        beta[4]=2.3038164e+00
        beta[5]=0
        beta[6]=5.0000000e-01
        beta[7]=0
            
        self.tm0_ig= self.compute_tm01_ig()
        self.IG= self.compute_ig_in()
        
        l0=np.squeeze(np.sqrt(9.81*0.33333*self.hm0)*self.tm0_ig)
        ksib=np.squeeze(self.sl2/(np.sqrt(self.IG/l0)))
        
        ksib=ksib-beta[5]
        ksib=np.maximum(ksib,0)
        ksibm=beta[4]*self.surfslope**beta[3]    # should increase with surfslope ie beta[4]>0
        cb=beta[2]*np.sqrt(self.surfslope)          # assume linear decrease ie beta[2]<0
        psibd=beta[0]*ksib**beta[1]
        psibr=(beta[0]*ksibm**beta[1] - 2.0)*(np.exp(-(ksib-ksibm)/(beta[2]))) + 2
        psib=psibd
        ind= np.argwhere(ksib>ksibm)
        psib[ind]=psibr[ind]
        
        v=self.IG*psib
        
        fac1= self.compute_dirspreadfac_lf(drspr)
        fac2= self.compute_directionfac_lf(phi)
        self.hm0_lf=v*fac1*fac2
        
    def compute_hm0_hf(self, drspr, phi):
        beta=[None]*8
        beta[0] =9.5635099e-01
        beta[1] =2.0143005e+00
        beta[2] =5.3602429e-01
        beta[3] =2.0000000e+00
        beta[4] =0.0000000e+00
        beta[5] =6.1856544e-01
        beta[6] =1.0000000e+00
        beta[7] =0.0000000e+00
       
        ksib=self.sl2/(np.sqrt(self.hm0/self.l1))
        
        self.hm0_hf=beta[0]*self.hm0*ksib**beta[1]*np.tanh((self.ksis+beta[4])**beta[5]/(beta[2]*ksib**beta[3]))
        
    def compute_tm01_ig(self):
        beta=[None]*5
        beta[0]=4.4021341e-07
        beta[1]=1.8635421e+00
        beta[2]=-4.2705433e-01
        beta[3]=7.2541023e-02
        beta[4]=2.0058478e+01
        
        tm0_ig=beta[0]+beta[1]*self.surfslope**beta[2]*self.steepness**beta[3]+beta[4]*self.surfslope
        return tm0_ig*self.tp  #self.tm0_ig

    def compute_ig_in(self):        
        beta=[None]*7
        beta[0]=2.2740842e+00
        beta[1]=1.0000000e+00
        beta[2]=5.0000000e-01
        beta[3]=2.7211454e+03
        beta[4]=2.0000000e+00
        beta[5]=1.7794945e+01
        beta[6]=1.8728433e-01

        return self.hm0*(beta[0]*self.surfslope**beta[2]*np.exp(-beta[5]*self.surfslope**beta[1]) 
                    + beta[6]*np.exp(-beta[3]*self.steepness**beta[4]) ) #self.IG

    def compute_directionfac_setup(self, phi):
        beta=[None]*5
        beta[0] = 1.4291 
        beta[1] = 0.0035124 
        beta[2] = 0.31891 
        beta[3] = 0.5 
        beta[4] = 0 
        
        return 1-beta[1]*self.steepness**beta[2]*np.absolute(phi)**beta[0] 

    def compute_dirspreadfac_setup(self, drspr):
        beta=[None]*5
        beta[0] = 0.031448 
        beta[1] = 0.69432 
        beta[2] = 0.66677 
        beta[3] = 0.5 
        beta[4] = 0 
        
        return np.exp(-beta[0]*drspr**beta[1] * (1.0 - np.tanh(beta[2]*self.ksis))) 
        
    def compute_dirspreadfac_lf(self, drspr):
        beta=[None]*5
        beta[0] = 0.047593
        beta[1] = 0.67228
        beta[2] = 0.50777
        beta[3] = 0.5
        beta[4] = 0
               
        return np.exp(-beta[0]*drspr**beta[1]* (1.0 - np.tanh(beta[2]*self.ksis)))
    
    def compute_directionfac_lf(self, phi):
        beta=[None]*5
        beta[0] = 0.40488
        beta[1] = 2.7073
        beta[2] = 0
        beta[3] = 0.5
        beta[4] = 0
        
        return (np.cos(phi*np.pi/180))**(beta[0]+beta[1]*self.surfslope)
        
    def compute_directionfac_hf(self, phi):
        beta=[None]*5
        beta[0] = 4.1355e-10
        beta[1] = 4.3559e-10
        beta[2] = 0 
        beta[3] = 0.5 
        beta[4] = 0 
        
        return (np.cos(phi*np.pi/180))**beta[0] 

    def compute_dirspreadfac_hf(self, drspr):
        beta=[None]*5
        beta[0] = 0.044544 
        beta[1] = 1 
        beta[2] = 21.1281 
        beta[3] = 0.5 
        beta[4] = 0 
        
        return np.exp(-beta[0]*drspr**beta[1] * (1.0 - np.tanh(beta[2]*self.ksis))) 
    
        
