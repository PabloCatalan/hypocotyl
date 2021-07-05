import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.ticker import FormatStrFormatter
prop = fm.FontProperties(fname='/usr/share/fonts/truetype/Mark Simonson - Proxima Nova Alt Regular-webfont.ttf')
import pylab
import numpy as np
import seaborn as sns
from scipy.integrate import odeint
import itertools
import subprocess
import pandas as pd

def get_key(df, key):
    kL=df[key].tolist()
    kM={}
    for k in kL:
        kM[k]=1
    kL=[]
    for k in kM:
        kL.append(k)
    return sorted(kL)

def read_data():
    DF=pd.read_csv('data/daylength_def.csv')
    mut=get_key(DF,'Genotype')
    Length=get_key(DF,'Daylength')
    Temp=get_key(DF,'Temperature')    
    avgdata={}
    stddata={}
    for m in mut:
        DFm=DF[DF.Genotype==m]
        m=m.replace(' ','_')
        avgdata[m]={}
        stddata[m]={}
        for d,t in itertools.product(Length,Temp):
            key=str(t)+'_'+str(d)
            DFnew=DFm[(DFm.Daylength==d) & (DFm.Temperature==t)]
            avgdata[m][key]=DFnew.Growth.mean()
            stddata[m][key]=DFnew.Growth.std()
    return avgdata, stddata

def is_day(t,Daylength):
    t1=t%24-Daylength
    return 1-np.heaviside(t1,1)


def elf3p(t,Daylength,pE1,pE2):  
    k0=5
    t2=t-Daylength
    t3=t-24.0
    
    if Daylength==0:
        return pE1+pE2
    elif Daylength==24:
        return pE1-pE2    
    else:
        SigT=2.0/(1.0+np.exp(-k0*t))
        SigT2=2.0/(1.0+np.exp(-k0*t2))
        SigT3=2.0/(1.0+np.exp(-k0*t3))
        return pE1-pE2*(-1.0+SigT-SigT2+SigT3)

def growth(y, t, Temp, Day, mut, pB28, kr22, kr28, pE122, pE128,
           pE222, pE228, dE, pPE22, pPE28, dP, kPC, dPB, pCL28,
           pCD, dC, pG, kG, pGP, pGE, pGB, pGH, pHC, mutBox,
           mutEox, mutPox, mutPko1, mutPko2, mutCox, mutCko1,
           mutCko2):
        
    #Variables
    B=y[0]#PHYB
    E=y[1]#ELF3
    P=y[2]#PIF
    C=y[3]#COP1
    G=y[4]#Hypocotyl
    #Parameters
    t1=t%24
    L=is_day(t1,Day)
    pB=10.0
    kr=kr22#0.232 datos de Casal
    pE1=pE122#adimensional
    pE2=pE222
    pP=1.0#adimensional
    pPE=pPE22
    pCL=1.0#adimensional
    mB=1.0#maximum PHYB value
    if Temp==28:
        pB=pB28
        kr=kr28#0.411 datos Casal
        pE1=pE128
        pE2=pE228
        pPE=pPE28
        pCL=pCL28
    if 'PHYBox' in mut:
        mB*=mutBox
    if 'ELF3ox' in mut:
        pE1*=mutEox
    if 'PIF4ox' in mut:
        pP*=mutPox
    if 'pif4' in mut:
        pP*=mutPko1
    if 'pifq' in mut:
        pP*=mutPko2
    if 'COP1' in mut:
        pCL*=mutCox
        pCD*=mutCox
    if 'cop1-4' in mut:
        pCL*=mutCko1
        pCD*=mutCko1
    if 'cop1-6' in mut:
        pCL*=mutCko2
        pCD*=mutCko2
    if 'hy5' in mut:
        pGH=0
    #Equations
    dBdt=pB*L*(mB-B)-kr*B
    dEdt=elf3p(t1,Day,pE1,pE2)-dE*E
    dPdt=pP/(1+pPE*E)-dP*P/(1+kPC*C)-dPB*P*B
    dCdt=pCL*L+pCD*(1-L)-dC*C
    dGdt=pG+kG*pGP*P/(1+pGP*P+pGE*E+pGB*B+pGH/(1+pHC*C))
    
    if 'elf3-8' in mut:
        dEdt=0
    if 'phyB' in mut:
        dBdt=0
    dydt=[dBdt, dEdt, dPdt, dCdt, dGdt]
    return dydt  
   
def model_results(Temp, Daylength, params, mutants):
    hypo={}
    tot={}
    for mut in mutants:
        hypo[mut]={}
        tot[mut]={}
        for T,D in itertools.product(Temp,Daylength):
                key=str(T)+'_'+str(D)
                time=np.linspace(0,120,500)
                y0=[0,0,0,0,0]
                sol=odeint(growth, y0, time, args=(T, D, mut, *params))
                R=sol[-1,4]
                hypo[mut][key]=R
                tot[mut][key]=sol
    return hypo, tot

def arrows(Temp, Daylength, params, mutants):
    S={}
    mut='Col'
    for T,D in itertools.product(Temp,Daylength):
        key=str(T)+'_'+str(D)
        time=np.linspace(0,120,500)
        y0=[0,0,0,0,0]
        sol=odeint(growth, y0, time, args=(T, D, mut, *params))
        S[key]=pd.DataFrame(sol, columns=['phyB', 'ELF3', 'PIF4', 'COP1', 'Growth'])
    return S

def decorrelate(Temp, Daylength, params, mutants):
    hypo={}
    tot={}
    for mut in mutants:
        hypo[mut]={}
        tot[mut]={}
        for T,D in itertools.product(Temp,Daylength):
                key=str(T)+'_'+str(D)
                time=np.linspace(0,120,500)
                y0=[0,0,0,0,0]
                sol=odeint(growth, y0, time, args=(T, D, mut, *params))
                R=sol[-1,4]
                hypo[mut][key]=R
                tot[mut][key]=sol
    return hypo, tot
