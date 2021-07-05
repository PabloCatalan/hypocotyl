#%%COMPUTE MODEL
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
from functions_paper import *
dir='/home/pcatalan/MEGA/plants_ode/'

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

#READ DARA
avgdata, stddata=read_data()

#SUFFIX FOR SIMULATIONS
suffix='01_paper'

#PARAMETERS
Temp=[22,28]
Daylength=[0,8,12,16,24]
Daylength2=np.linspace(0,24,49)#range(0,25)
Daylength2=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4]+list(range(4,25,2))
#Daylength2=list(Daylength)
params=[]
with open(dir+'results/parameters_'+suffix+'.txt', 'r') as f:
    for line in f:
        w=line.split()
        params.append(float(w[0]))       
        
#%%FIGURE 1
import time

mutants=['Col', 'pif4', 'pifq', 'PIF4ox',
         'phyB', 'PHYBox', 'elf3-8', 'ELF3ox',
         'COP1ox', 'cop1-4', 'ELF3ox_cop1-4', 'elf3-8_cop1-4',
         'phyB_cop1-4', 'phyB_COP1ox', 'phyB_elf3-8', 'PHYBox_COP1ox']
# mutants=['Col', 'phyB', 'cop1-6']
# mutants=['Col', 'PHYBox', 'ELF3ox', 'PHYBox_ELF3ox']
#SIMULATION IN C++
#for mut in mutants:
#    com='./new_model_paper '+str(mut)+' '+suffix
#    subprocess.call([com], shell=True)
##READ PREDICTIONS FROM C++ MODEL
#hypo={}
#for mut in mutants:
#    hypo[mut]={}
#    with open(dir+'results/nmpc++_results_'+mut+'_'+suffix+'.txt','r') as f:
#        for line in f:
#            w=line.split()
#            T=w[0]
#            D=w[1]
#            g=float(w[2])
#            key=T+'_'+str(D)
#            hypo[mut][key]=g
#SIMULATION IN PYTHON
hypo_python, tot_python=model_results(Temp,Daylength2,params,mutants)

#%%
mutants=['Col', 'pif4', 'pifq', 'PIF4ox',
         'phyB', 'PHYBox', 'elf3-8', 'ELF3ox',
         'COP1ox', 'cop1-4', 'ELF3ox_cop1-4', 'elf3-8_cop1-4']
#WRITE DATA
with pd.ExcelWriter(dir+'results/fig1_data.xls') as writer:
    for mut in mutants:
        W=[]
        for (D,T) in itertools.product(Daylength2,Temp):
            key=str(T)+'_'+str(D)
            try:
                W.append([D,T,avgdata[mut][key],stddata[mut][key],hypo_python[mut][key]])
            except:
                W.append([D,T,'NaN','NaN',hypo_python[mut][key]])
        DFdata=pd.DataFrame.from_records(W, columns=['Daylength', 'Temperature', 'Average Growth (experimental)', 'Standard Deviation Growth (experimental)', 'Model prediction'])
        DFdata.to_excel(writer, sheet_name=mut, float_format='%.3f',index=False)

#%%PREDICTIONS
pred=['phyB_cop1-4', 'phyB_COP1ox', 'phyB_elf3-8', 'PHYBox_COP1ox']
#WRITE DATA
with pd.ExcelWriter(dir+'results/predicted_data.xls') as writer:
    for mut in pred:
        W=[]
        for (D,T) in itertools.product(Daylength2,Temp):
            key=str(T)+'_'+str(D)
            W.append([D,T,'NaN','NaN',hypo_python[mut][key]])
        DFdata=pd.DataFrame.from_records(W, columns=['Daylength', 'Temperature', 'Average Growth (experimental)', 'Standard Deviation Growth (experimental)', 'Model prediction'])
        DFdata.to_excel(writer, sheet_name=mut, float_format='%.3f',index=False)

# for D in Daylength2:
#         key22='22_'+str(D)
#         key28='28_'+str(D)
#         hp22.append(hypo_python[mut][key22])
#         hp28.append(hypo_python[mut][key28])  

#%%
#PLOT
fig=plt.figure(figsize=(7,6))
ncols=4
nrows=4
for i1,mut in enumerate(mutants):
    ax=fig.add_subplot(nrows,ncols,i1+1)
    #PLOT SIMULATIONS C++
#    h22=[]
#    h28=[]
#    for D in Daylength2:
#        key22='22_'+str(D)
#        key28='28_'+str(D)
#        h22.append(hypo[mut][key22])
#        h28.append(hypo[mut][key28])
#    ax.plot(Daylength2, h22, '--k')
#    ax.plot(Daylength2, h28, '--r')
    #PLOT SIMULATIONS PYTHON
    hp22=[]
    hp28=[]
    for D in Daylength2:
        key22='22_'+str(D)
        key28='28_'+str(D)
        hp22.append(hypo_python[mut][key22])
        hp28.append(hypo_python[mut][key28])  
    ax.plot(Daylength2, hp22, 'k', label='22 ºC')
    ax.plot(Daylength2, hp28, 'r', label='28 ºC')
    #DATA
    d22=[]
    d28=[]
    s22=[]
    s28=[]
    for D in Daylength:
        key22='22_'+str(D)
        key28='28_'+str(D)
        if mut in avgdata:
            d22.append(avgdata[mut][key22])
            d28.append(avgdata[mut][key28])
            s22.append(stddata[mut][key22])
            s28.append(stddata[mut][key28])
    if mut in avgdata:
        ax.errorbar(Daylength, d22, yerr=s22, fmt='o', color='k')
        ax.errorbar(Daylength, d28, yerr=s28, fmt='o', color='r')
    if mut=='Col':
        ax.set_title(mut, size=10)
    else:
        ax.set_title(mut, style='italic', size=10)
    ax.set_ylim([0,20])
    if i1==0:
        ax.legend(loc='upper right', frameon=False)
    if i1>7:
        ax.set_xlabel('Daylength (hours)', size=10)
    if i1%4==0:
        ax.set_ylabel('growth (mm)', size=10)
    ax.set_xticks([0,4,8,12,16,20,24])
    ax.set_xticklabels([0,4,8,12,16,20,24], size=5)
    ax.set_yticklabels(ax.get_yticks(), size=5)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
fig.tight_layout()
fig.savefig('figures/fig1_'+suffix+'.pdf', bbox_inches='tight')

#%%SUPP FIGURE 1 (EXTRA MUTANTS)
extraM=['cop1-6', 'phyB_cop1-6', 'hy5']
# extraM=['Col', 'PHYBox', 'ELF3ox', 'PHYBox_ELF3ox']
hypo_supp, tot_supp=model_results(Temp,Daylength2,params,extraM)
fig=plt.figure(figsize=(8,2))
ncols=len(extraM)
nrows=1
for i1,mut in enumerate(extraM):
    ax=fig.add_subplot(nrows,ncols,i1+1)
    #PLOT SIMULATIONS PYTHON
    hp22=[]
    hp28=[]
    for D in Daylength2:
        key22='22_'+str(D)
        key28='28_'+str(D)
        hp22.append(hypo_supp[mut][key22])
        hp28.append(hypo_supp[mut][key28])  
    ax.plot(Daylength2, hp22, 'k', label='22 ºC')
    ax.plot(Daylength2, hp28, 'r', label='28 ºC')
    #DATA
    d22=[]
    d28=[]
    s22=[]
    s28=[]
    for D in Daylength:
        key22='22_'+str(D)
        key28='28_'+str(D)
        if mut in avgdata:
            d22.append(avgdata[mut][key22])
            d28.append(avgdata[mut][key28])
            s22.append(stddata[mut][key22])
            s28.append(stddata[mut][key28])
    if mut in avgdata:
        ax.errorbar(Daylength, d22, yerr=s22, fmt='o', color='k')
        ax.errorbar(Daylength, d28, yerr=s28, fmt='o', color='r')
    ax.set_title(mut, style='italic', size=10)
    ax.set_ylim([0,20])
    if i1==0:
        ax.legend(loc='upper right', prop=prop, frameon=False)
    if i1>=0:
        ax.set_xlabel('Daylength (hours)', fontproperties=prop, size=10)
    if i1%4==0:
        ax.set_ylabel('growth (mm)', fontproperties=prop, size=10)
    ax.set_xticks([0,4,8,12,16,20,24])
    ax.set_xticklabels([0,4,8,12,16,20,24], size=10, fontproperties=prop)
    ax.set_yticklabels(ax.get_yticks(), size=10, fontproperties=prop)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
fig.tight_layout()
fig.savefig('figures/suppfig1'+suffix+'.pdf', bbox_inches='tight')
#%%WRITE DATA SUPPFIG1
extraM=['cop1-6', 'phyB_cop1-6', 'hy5']
#WRITE DATA
with pd.ExcelWriter(dir+'results/suppfig1_data.xls') as writer:
    for mut in extraM:
        W=[]
        for (D,T) in itertools.product(Daylength2,Temp):
            key=str(T)+'_'+str(D)
            try:
                W.append([D,T,avgdata[mut][key],stddata[mut][key],hypo_supp[mut][key]])
            except:
                W.append([D,T,'NaN','NaN',hypo_supp[mut][key]])
        DFdata=pd.DataFrame.from_records(W, columns=['Daylength', 'Temperature', 'Average Growth (experimental)', 'Standard Deviation Growth (experimental)', 'Model prediction'])
        DFdata.to_excel(writer, sheet_name=mut, float_format='%.3f',index=False)


#     DFd=DF[DF['D']==24]
#     DFd.to_excel(writer, sheet_name='Light', float_format='%.3f',index=False)
# for mut in extraM:
#     for (D,T) in itertools.product(Daylength,Temp):
#         key=str(T)+'_'+str(D)
#         W.append([mut,D,T,avgdata[mut][key],stddata[mut][key],hypo_supp[mut][key]])
# DFdata=pd.DataFrame.from_records(W, columns=['Genotype', 'Daylength', 'Temperature', 'Average Growth (experimental)', 'Standard Deviation Growth (experimental)', 'Model prediction'])
# DFdata.to_excel(dir+'results/suppfig1_data.xls',float_format='%.3f')
#%%FIG 2 (PROTEIN DYNAMICS)
time=np.linspace(0,120,500)
mut='Col'
B22=tot_python[mut]['22_8'][:,0]
B28=tot_python[mut]['28_8'][:,0]
E22=tot_python[mut]['22_8'][:,1]
E28=tot_python[mut]['28_8'][:,1]
P22=tot_python[mut]['22_8'][:,2]
P28=tot_python[mut]['28_8'][:,2]
C22=tot_python[mut]['22_8'][:,3]
C28=tot_python[mut]['28_8'][:,3]
G22=tot_python[mut]['22_8'][:,4]
G28=tot_python[mut]['28_8'][:,4]

DF=pd.DataFrame({'Time':time, 'phyb22':B22, 'phyb28':B28,
                 'ELF322':E22, 'ELF328':E28, 'PIF422':P22,
                 'PIF428':P28, 'COP122':C22, 'COP128':C28})
DF.to_csv(dir+'results/short_day_proteins.csv', float_format='%.3f', index=False)
DF.to_excel(dir+'results/short_day_proteins.xls', float_format='%.3f', index=False)

def modsavefigs(time, b22, b28, Daylength, Title):
    xmin=24
    xmax=120
    mask=(time>=xmin) & (time<=xmax)
    Time=time[mask]
    B22=b22[mask]
    B28=b28[mask]
    fig=plt.figure()
    ax=plt.gca()
    ax.plot(Time,B22,'k', label='22 ºC')
    ax.plot(Time,B28,'r', label='28 ºC')
    for days in range(0,5):
        time_day=np.linspace(days*24+Daylength,(days*24)+24, 100)
        ax.fill_between(time_day, 0, 10000, facecolor='grey', alpha=0.5)
    ax.set_ylim([0.9*min(min(B22),min(B28)), 1.10*max(max(B22),max(B28))])
    ax.set_xlabel('ZT (h)', fontproperties=prop, size=15)
    ax.set_ylabel(Title, fontproperties=prop, size=15)
    ax.set_xticklabels(ax.get_xticks(), size=10, fontproperties=prop)
    ax.set_yticklabels(ax.get_yticks(), size=10, fontproperties=prop)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.legend(loc='upper right',frameon=False, prop=prop)
    ax.set_xlim([xmin, xmax])
    if Title=='COP1':
        ax.set_yscale('log')
    fig.savefig('figures/fig2_'+Title+'_'+suffix+'.pdf', bbox_inches='tight')
modsavefigs(time, B22, B28, 8, 'PHYB')
modsavefigs(time, E22, E28, 8, 'ELF3')
modsavefigs(time, P22, P28, 8, 'PIF4')
modsavefigs(time, C22, C28, 8, 'COP1')
modsavefigs(time, G22, G28, 8, 'growth (mm)')

#%%FIG 2 (PROTEIN DYNAMICS In long day)
time=np.linspace(0,120,500)
B22=tot_python['Col']['22_16'][:,0]
B28=tot_python['Col']['28_16'][:,0]
E22=tot_python['Col']['22_16'][:,1]
E28=tot_python['Col']['28_16'][:,1]
P22=tot_python['Col']['22_16'][:,2]
P28=tot_python['Col']['28_16'][:,2]
C22=tot_python['Col']['22_16'][:,3]
C28=tot_python['Col']['28_16'][:,3]
G22=tot_python['Col']['22_16'][:,4]
G28=tot_python['Col']['28_16'][:,4]

DF=pd.DataFrame({'Time':time, 'phyb22':B22, 'phyb28':B28,
                 'ELF322':E22, 'ELF328':E28, 'PIF422':P22,
                 'PIF428':P28, 'COP122':C22, 'COP128':C28})
DF.to_csv(dir+'results/long_day_proteins.csv', float_format='%.3f', index=False)
DF.to_excel(dir+'results/long_day_proteins.xls', float_format='%.3f', index=False)

def modsavefigs(time, b22, b28, Daylength, Title):
    xmin=24
    xmax=120
    mask=(time>=xmin) & (time<=xmax)
    Time=time[mask]
    B22=b22[mask]
    B28=b28[mask]
    fig=plt.figure()
    ax=plt.gca()
    ax.plot(Time,B22,'k', label='22 ºC')
    ax.plot(Time,B28,'r', label='28 ºC')
    for days in range(0,5):
        time_day=np.linspace(days*24+Daylength,(days*24)+24, 100)
        ax.fill_between(time_day, 0, 10000, facecolor='grey', alpha=0.5)
    ax.set_ylim([0.9*min(min(B22),min(B28)), 1.10*max(max(B22),max(B28))])
    ax.set_xlabel('ZT (h)', fontproperties=prop, size=15)
    ax.set_ylabel(Title, fontproperties=prop, size=15)
    ax.set_xticklabels(ax.get_xticks(), size=10, fontproperties=prop)
    ax.set_yticklabels(ax.get_yticks(), size=10, fontproperties=prop)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.legend(loc='upper right',frameon=False, prop=prop)
    ax.set_xlim([xmin, xmax])
    fig.savefig('figures/SUPPfig2_'+Title+'_'+suffix+'.pdf', bbox_inches='tight')
modsavefigs(time, B22, B28, 16, 'PHYB')
modsavefigs(time, E22, E28, 16, 'ELF3')
modsavefigs(time, P22, P28, 16, 'PIF4')
modsavefigs(time, C22, C28, 16, 'COP1')
modsavefigs(time, G22, G28, 16, 'growth (mm)')
#%%FIG 3 (EFFECT OF COP1)
DF=pd.read_csv(dir+'results/effect_of_cop1_'+suffix+'.csv')
with pd.ExcelWriter(dir+'results/fig3_data.xls') as writer:
    DFd=DF[DF['D']==0]
    DFd.to_excel(writer, sheet_name='Dark', index=False)
    DFd=DF[DF['D']==24]
    DFd.to_excel(writer, sheet_name='Light', index=False)
# DF.to_excel(dir+'results/fig3_data.xls', float_format='%.5f')
#%%
fig=plt.figure()
#DARK
ax=fig.add_subplot(1,2,1)
DFd=DF[DF['D']==0]
TC=['k', 'r']
COP1m=['COP1ox', 'cop1-4', 'Col', 'hy5']
COP1ex=[params[28], params[29],  1.0, 2000]
SL=['^', 's', 'o', 'x', '^']
for t1,t in enumerate([22,28]):
    DFc=DFd[DFd['T']==t]
    ax.plot(DFc['COP1'], DFc['Growth'], color=TC[t1], label=str(t)+' ºC')
    ax.set_xscale('log')
    ax.set_xlim([0.01,3000])
    ax.set_ylim([0,16])
    ax.set_xlabel('relative COP1 expression', fontproperties=prop, size=15)
    ax.set_ylabel('growth (mm)', fontproperties=prop, size=15)
    ax.set_xticklabels(ax.get_xticks(), size=10, fontproperties=prop)
    ax.set_yticklabels(ax.get_yticks(), size=10, fontproperties=prop)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.set_title('Dark', fontproperties=prop, fontsize=15)
    ax.legend(loc='lower right', prop=prop, frameon=False)
    #PLOT MUTANTS
    for c1,c in enumerate(COP1m):
        sT=str(t)
        y=avgdata[COP1m[c1]][sT+'_0']
        yerr=stddata[COP1m[c1]][sT+'_0']
        ax.errorbar(COP1ex[c1],y,yerr=yerr,fmt=SL[c1],color='w', markersize=8,ecolor=TC[t1],markeredgecolor=TC[t1])
        print(0,t,COP1ex[c1],y,yerr)
#    xp=[params[28], params[29], params[30], 1.0]
#    sT=str(t)
#    yp=[avgdata['COP1ox'][sT+'_0'],avgdata['cop1-4'][sT+'_0'],avgdata['cop1-6'][sT+'_0'],avgdata['Col'][sT+'_0']]
#    yerr=[stddata['COP1ox'][sT+'_0'],stddata['cop1-4'][sT+'_0'],stddata['cop1-6'][sT+'_0'],stddata['Col'][sT+'_0']]
#    ax.errorbar(xp,yp,yerr=yerr, fmt='o', color=TC[t1])
#LIGHT
ax=fig.add_subplot(1,2,2)
DFd=DF[DF['D']==24]
TC=['k', 'r']
for t1,t in enumerate([22,28]):
    DFc=DFd[DFd['T']==t]
    ax.plot(DFc['COP1'], DFc['Growth'], color=TC[t1], label='_nolegend_')
    ax.set_xscale('log')
    ax.set_xlim([0.01,3000])
    ax.set_ylim([0,16])
    ax.set_xlabel('relative COP1 expression', fontproperties=prop, size=15)
    ax.set_ylabel('growth (mm)', fontproperties=prop, size=15)
    ax.set_xticklabels(ax.get_xticks(), size=10, fontproperties=prop)
    ax.set_yticklabels(ax.get_yticks(), size=10, fontproperties=prop)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.set_title('CWL', fontproperties=prop, fontsize=15)
    #PLOT MUTANTS
    for c1,c in enumerate(COP1m):
        sT=str(t)
        y=avgdata[COP1m[c1]][sT+'_24']
        yerr=stddata[COP1m[c1]][sT+'_24']
        if t==22:
            ax.errorbar(COP1ex[c1],y,yerr=yerr,fmt=SL[c1],color='w',markersize=8,ecolor=TC[t1],markeredgecolor=TC[t1],label=c)
            print(24,t,COP1ex[c1],y,yerr)
        else:
            ax.errorbar(COP1ex[c1],y,yerr=yerr,fmt=SL[c1],color='w',markersize=8,ecolor=TC[t1],markeredgecolor=TC[t1])
            print(24,t,COP1ex[c1],y,yerr)
    ax.legend(prop=prop, frameon=False)
#    xp=[params[28], params[29], params[30], 1.0]
#    sT=str(t)
#    yp=[avgdata['COP1ox'][sT+'_24'],avgdata['cop1-4'][sT+'_24'],avgdata['cop1-6'][sT+'_24'],avgdata['Col'][sT+'_24']]
#    yerr=[stddata['COP1ox'][sT+'_24'],stddata['cop1-4'][sT+'_24'],stddata['cop1-6'][sT+'_24'],stddata['Col'][sT+'_24']]
#    ax.errorbar(xp,yp,yerr=yerr, fmt='o', color=TC[t1])
fig.tight_layout()
fig.savefig(dir+'figures/fig3_cop1_expression.pdf', bbox_inches='tight')

#%%SUPPFIG DATA FOR FITTING
DF=pd.read_csv(dir+'data/daylength_def.csv')
mutants=sorted(get_key(DF,'Genotype'),key=str.casefold)
mutants.remove('lux')
mutants.remove('elf4')
mutants.remove('det1-1')
fig=plt.figure(figsize=(16,12))
maxcols=5
maxrows=int(np.ceil(len(mutants)/maxcols))
for i,mut in enumerate(mutants):
    ax=fig.add_subplot(maxrows,maxcols,i+1)
    M=DF[DF.Genotype==mut]
    TC=['k','r']
    for t1,T in enumerate([22,28]):
        MT=M[M.Temperature==T]
        D=MT.Daylength+np.random.normal(0,0.5, len(MT))
        G=MT.Growth
        ax.scatter(D,G,color=TC[t1],label=str(T)+' ºC')
    from collections import OrderedDict
    handles, labels = ax.get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    if i==0:
        ax.legend(by_label.values(), by_label.keys(),loc='upper right', prop=prop, frameon=False)
    ax.set_title(mut, style='italic', size=15)
    ax.set_ylim([0,20])
    if i>=10:
        ax.set_xlabel('Daylength (hours)', fontproperties=prop, size=15)
    else:
        ax.set_xlabel('')
    if i%5==0:
        ax.set_ylabel('growth (mm)', fontproperties=prop, size=15)
    else:
        ax.set_ylabel('')
    ax.set_xticks([0,4,8,12,16,20,24])
    ax.set_xticklabels(ax.get_xticks(), size=10, fontproperties=prop)
    ax.set_yticklabels(ax.get_yticks(), size=10, fontproperties=prop)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
fig.tight_layout()
fig.savefig('figures/suppfig_data_growth_'+suffix+'.pdf', bbox_inches='tight')

#%%
EX=pd.read_csv(dir+'data/expression_def.csv', sep='\t')
T=np.arange(0,24,1)
T2=np.arange(0,24,0.001)
fig=plt.figure(figsize=(6,3))
#COL DATA
ax=fig.add_subplot(1,2,1)
ax.plot(T,EX['Col 22'], 'o-k', label='22 ºC', markersize=4)
ax.plot(T,EX['Col 28'], 'o-r', label='28 ºC', markersize=4)
D=np.heaviside(T2-8,0)
ax.fill_between(T2,0*D,10*D,facecolor='grey', alpha=0.5)
ax.set_title('WT', size=15)
ax.set_xticks([0,4,8,12,16,20,24])
ax.set_xticklabels(ax.get_xticks(), size=10, fontproperties=prop)
ax.set_yticklabels(ax.get_yticks(), size=10, fontproperties=prop)
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.set_xlabel('ZT(h)', fontproperties=prop, size=15)
ax.set_ylabel('Relative luminiscence', fontproperties=prop, size=15)
#cop1-6 DATA
# ax=fig.add_subplot(1,3,2)
# ax.plot(T,EX['cop1-6 22'], '-k', label='22 ºC')
# ax.plot(T,EX['cop1-6 28'], '-r', label='28 ºC')
# D=np.heaviside(T-8,0)
# ax.fill_between(T,0*D,10*D,facecolor='grey', alpha=0.5)
# ax.set_title('cop1-6', style='italic', size=15)
# ax.set_xticks([0,4,8,12,16,20,24])
# ax.set_xticklabels(ax.get_xticks(), size=10, fontproperties=prop)
# ax.set_yticklabels(ax.get_yticks(), size=10, fontproperties=prop)
# ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
# ax.set_xlabel('ZT(h)', fontproperties=prop, size=15)
# ax.set_ylabel('Relative luminiscence', fontproperties=prop, size=15)
#phyB DATA
ax=fig.add_subplot(1,2,2)
ax.plot(T,EX['phyB 22'], 'o-k', label='22 ºC', markersize=4)
ax.plot(T,EX['phyB 28'], 'o-r', label='28 ºC', markersize=4)
D=np.heaviside(T2-8,0)
ax.fill_between(T2,0*D,10*D,facecolor='grey', alpha=0.5)
ax.set_title('phyB-9', style='italic', size=15)
ax.set_xticks([0,4,8,12,16,20,24])
ax.set_xticklabels(ax.get_xticks(), size=10, fontproperties=prop)
ax.set_yticklabels(ax.get_yticks(), size=10, fontproperties=prop)
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.set_xlabel('ZT(h)', fontproperties=prop, size=15)
ax.set_ylabel('Relative luminiscence', fontproperties=prop, size=15)

fig.tight_layout()
fig.savefig('figures/suppfig_data_expression_'+suffix+'.pdf', bbox_inches='tight')

#%%SUPPFIG DATA FOR FITTING VIOLINPLOT
DF=pd.read_csv(dir+'data/daylength_def.csv')
mutants=sorted(get_key(DF,'Genotype'),key=str.casefold)
mutants.remove('lux')
mutants.remove('elf4')
mutants.remove('det1-1')
# mutants=['Col', 'pif4', 'pifq', 'PIF4ox',
#          'phyB', 'PHYBox', 'elf3-8', 'ELF3ox',
#          'COP1ox', 'cop1-4', 'ELF3ox cop1-4', 'elf3-8 cop1-4']
fig=plt.figure(figsize=(16,12))
maxcols=5
maxrows=int(np.ceil(len(mutants)/maxcols))
for i,mut in enumerate(mutants):
    ax=fig.add_subplot(maxrows,maxcols,i+1)
    M=DF[DF.Genotype==mut]
    TC=['k','r']
    TP=['binary', 'Reds']
#    for t1,T in enumerate([22,28]):
#        MT=M[M.Temperature==T]
    D=M.Daylength
    G=M.Growth
    ax1=sns.violinplot(x='Daylength', y='Growth', hue='Temperature', 
                       data=M, palette=sns.color_palette(['darkgrey','r'], desat=0.5),
                       )#, linewidth=1)
    if i>0:
        ax1.get_legend().remove()#ax1.legend()
    
    # #PLOT SIMULATIONS PYTHON
    # hp22=[]
    # hp28=[]
    # for D in Daylength2:
    #     newmut=mut.replace(' ','_')
    #     key22='22_'+str(D)
    #     key28='28_'+str(D)
    #     hp22.append(hypo_python[newmut][key22])
    #     hp28.append(hypo_python[newmut][key28])
    # D3=np.array(Daylength2)/6
    # ax.plot(D3-0.25, hp22, 'grey', label='22 ºC')
    # ax.plot(D3+0.25, hp28, 'r', label='28 ºC')
#    from collections import OrderedDict
#    handles, labels = ax.get_legend_handles_labels()
#    by_label = OrderedDict(zip(labels, handles))
#    if i==0:
#        ax.legend(by_label.values(), by_label.keys(),loc='upper right', prop=prop, frameon=False)
    if mut=='Col':
        ax.set_title('WT', size=15)
    else:
        ax.set_title(mut, style='italic', size=15)
    ax.set_ylim([0,20])
    if i>=10:
        ax.set_xlabel('Daylength (hours)', fontproperties=prop, size=15)
    else:
        ax.set_xlabel('')
    if i%5==0:
        ax.set_ylabel('growth (mm)', fontproperties=prop, size=15)
    else:
        ax.set_ylabel('')
#    ax.set_xticks([0,4,8,12,16,20,24])
#    ax.set_xticklabels(ax.get_xticks(), size=10, fontproperties=prop)
#    ax.set_yticklabels(ax.get_yticks(), size=10, fontproperties=prop)
#    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
fig.tight_layout()
fig.savefig('figures/suppfig_data_growth_'+suffix+'_violin.pdf', bbox_inches='tight')

#%%SUPPFIG DATA FOR FITTING VIOLINPLOT
DF=pd.read_csv(dir+'data/daylength_def.csv')
mutants=sorted(get_key(DF,'Genotype'),key=str.casefold)
mutants.remove('lux')
mutants.remove('elf4')
mutants.remove('det1-1')
fig=plt.figure(figsize=(16,12))
maxcols=5
maxrows=int(np.ceil(len(mutants)/maxcols))
mutants=['PIF4ox', 'phyB']
for i,mut in enumerate(mutants):
    ax=fig.add_subplot(maxrows,maxcols,i+1)
    M=DF[DF.Genotype==mut]
    TC=['k','r']
    TP=['Greys', 'Reds']
#    for t1,T in enumerate([22,28]):
#        MT=M[M.Temperature==T]
    D=M.Daylength
    G=M.Growth
    ax1=sns.violinplot(x='Daylength', y='Growth', hue='Temperature', 
                       data=M, palette=sns.color_palette(['grey','r'], desat=0.5),
                       scale='count')#, linewidth=1)
    
    #PLOT SIMULATIONS PYTHON
    hp22=[]
    hp28=[]
    for D in Daylength2:
        key22='22_'+str(D)
        key28='28_'+str(D)
        hp22.append(hypo_python[mut][key22])
        hp28.append(hypo_python[mut][key28])  
    ax.plot(Daylength2, hp22, 'k', label='22 ºC')
    ax.plot(Daylength2, hp28, 'r', label='28 ºC')
#    from collections import OrderedDict
#    handles, labels = ax.get_legend_handles_labels()
#    by_label = OrderedDict(zip(labels, handles))
#    if i==0:
#        ax.legend(by_label.values(), by_label.keys(),loc='upper right', prop=prop, frameon=False)
    if mut=='Col':
        ax.set_title(mut, size=15)
    else:
        ax.set_title(mut, style='italic', size=15)
    ax.set_ylim([0,20])
    if i>=10:
        ax.set_xlabel('Daylength (hours)', fontproperties=prop, size=15)
    else:
        ax.set_xlabel('')
    if i%5==0:
        ax.set_ylabel('growth (mm)', fontproperties=prop, size=15)
    else:
        ax.set_ylabel('')
#    ax.set_xticks([0,4,8,12,16,20,24])
#    ax.set_xticklabels(ax.get_xticks(), size=10, fontproperties=prop)
#    ax.set_yticklabels(ax.get_yticks(), size=10, fontproperties=prop)
#    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
fig.tight_layout()
fig.savefig('figures/suppfig_data_growth_'+suffix+'_violin_'+str(mutants)+'.pdf', bbox_inches='tight')
