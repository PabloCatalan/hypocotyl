#%%
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
prop = fm.FontProperties(fname='/usr/share/fonts/truetype/Mark Simonson - Proxima Nova Alt Regular-webfont.ttf')
fname='/usr/share/fonts/truetype/Mark Simonson - Proxima Nova Alt Regular-webfont.ttf'
from matplotlib.ticker import FormatStrFormatter
import pylab
import numpy as np
import seaborn as sns
from scipy.integrate import odeint
import itertools
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pandas as pd
from operator import itemgetter
from scipy import stats

#GET KEY FROM DATAFRAME
def get_key(df, key):
    kL=df[key].tolist()
    kM={}
    for k in kL:
        kM[k]=1
    kL=[]
    for k in kM:
        kL.append(k)
    return sorted(kL)

def plot_heatmap(X1,Y1,R,mut,ylabel,vmax,X,Y):
    text=mut.split()
    prot=text[0][:-2]
    fondo=text[1]
    if len(text)>2:
        fondo=fondo+' '+text[2]
    fig=plt.figure()
    ax=plt.gca()
    cmap='coolwarm'
    if 'growth' in ylabel:
        cmap='summer_r'#'YlGn'
    I=ax.pcolormesh(X,Y,R, vmax=vmax, vmin=0, rasterized=True, cmap=cmap)#'RdBu_r')
    ax.set_xscale('log')
    # ax.set_xlim([10**(-2.5),10**(2.5)])
    cr=0.2
    C=np.arange(0.1,vmax,cr)
    if 'growth' in ylabel: 
        cr=1
        C=np.arange(0,vmax+1,cr)
    CS=ax.contour(X,Y,R,C, linewidths=0.3, colors='k')#cmap='binary')
    CL=ax.clabel(CS, C, inline=True, fmt='%.1f', fontsize=7, use_clabeltext=True)
    cbar=fig.colorbar(I, ax=ax)
    cbar.ax.set_yticklabels(cbar.ax.get_yticks(), fontname='Helvetica')
    cbar.ax.set_ylabel(ylabel, fontname='Helvetica', fontsize=15)
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_ylabel('Daylength', fontname='Helvetica', fontsize=15)
    ax.set_xlabel(prot+'  activity', fontname='Helvetica', fontsize=15)
    ax.set_yticks(range(0,25,4))
    ax.set_yticklabels(ax.get_yticks(), fontname='Helvetica')
    ax.set_xticklabels(ax.get_xticks(), fontname='Helvetica')
    ax.set_title(fondo)
    # ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    plt.minorticks_off()
    ax.set_ylim([0,24])
    fig.savefig(dir+'figures/effect_of_'+prot+'_sobre_fondo_'+fondo+'_'+ylabel+'_heatmap.pdf', bbox_inches='tight')
    fig.savefig(dir+'figures/effect_of_'+prot+'_sobre_fondo_'+fondo+'_'+ylabel+'_heatmap.png', bbox_inches='tight')
    plt.close()

dir='/home/pcatalan/MEGA/plants_ode/'
M=['COP1ox WT', 'COP1ox phyb', 'COP1ox elf3-8', 'COP1ox phyb elf3-8','COP1ox pifq',
	'PHYBox WT', 'PHYBox cop1-4', 'PHYBox elf3-8', 'PHYBox cop1-4 elf3-8', 'PHYBox pifq',
	'ELF3ox WT', 'ELF3ox phyb', 'ELF3ox cop1-4', 'ELF3ox cop1-4 phyb', 'ELF3ox pifq',
    'PIF4ox WT', 'PIF4ox phyb', 'PIF4ox cop1-4', 'PIF4ox elf3-8']
# M=['COP1ox pifq', 'PHYBox pifq','ELF3ox pifq',
   # 'PIF4ox WT', 'PIF4ox phyb', 'PIF4ox cop1-4', 'PIF4ox elf3-8']
# M=['PHYBox WT']
for mut in M:
    print(mut)
    DF=pd.read_csv(dir+'results/heatmap_'+mut+'.csv')
    X=get_key(DF,'Prot')    
    Y=get_key(DF,'D')
    X=np.array(X)
    R=np.zeros((len(Y),len(X)))
    R1=np.zeros((len(Y),len(X)))
    R2=np.zeros((len(Y),len(X)))
    for i,x in enumerate(X):
        for j,y in enumerate(Y):
            A=DF[(DF['Prot']==x) & (DF['D']==y)]
            A22=A[A['T']==22]
            A28=A[A['T']==28]
            R[j,i]=A28['Growth'].mean()-A22['Growth'].mean()
            R1[j,i]=A28['Growth'].mean()
            R2[j,i]=A22['Growth'].mean()
    # np.savetxt(dir+'results/phyB_matrix.txt', R, fmt='%.3f')
    # np.savetxt(dir+'results/phyB_xaxis.txt', X, fmt='%.6f')
    # np.savetxt(dir+'results/phyB_yaxis.txt', Y, fmt='%.1f')
    yshift=0.5#0.25
    yadd=24+yshift
    Y1=np.array(Y)-yshift
    Y1=list(Y)+[yadd]
    xshift=0.1#0.05
    xadd=10**(3+xshift)
    X1=np.log10(X)-xshift
    X1=list(10**X1)+[xadd]
    # X=np.log10(X)
    
    vmax=max(R1.max(),R2.max())
    vmax=14
    Vmax=7.0
    plot_heatmap(X1,Y1,R,mut,'Hypocotyl Thermoelongation',Vmax,X,Y)
    plot_heatmap(X1,Y1,R1,mut,'growth at 28',vmax,X,Y)
    plot_heatmap(X1,Y1,R2,mut,'growth at 22',vmax,X,Y)
    
    #CURVES
    D=[0,4,8,12,16,20,24]
    vals=np.linspace(0,1,len(D))#each class will have one color
    cmap=plt.cm.colors.ListedColormap(plt.cm.viridis(vals))
    fig=plt.figure()
    ax=plt.gca()
    text=mut.split()
    prot=text[0][:-2]
    fondo=text[1]
    if len(text)>2:
        fondo=fondo+' '+text[2]
    for i,d in enumerate(D):
        A=DF[DF['D']==d]
        A22=A[A['T']==22]
        A28=A[A['T']==28]
        x=np.array(A22['Prot'])
        y=np.array(A28['Growth'])-np.array(A22['Growth'])
        DFnew=pd.DataFrame({'Prot':x, 'Thermo':y})
        DFnew.to_excel(dir+'results/effect_of_'+prot+'_sobre_fondo_'+fondo+'_'+str(d)+'.xls', float_format='%.6f', index=False)
        ax.plot(x,y,color=cmap(i/len(D)), label='D='+str(d))
    ax.legend()
    ax.set_ylabel('Thermomorphogenetic response', fontproperties=prop, fontsize=15)
    ax.set_xlabel('phyB activity', fontproperties=prop, fontsize=15)
    ax.set_xscale('log')
    fig.savefig(dir+'figures/effect_of_'+prot+'_sobre_fondo_'+fondo+'_curves.pdf', bbox_inches='tight')

#%%RELATIVE HEATMAPS
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
prop = fm.FontProperties(fname='/usr/share/fonts/truetype/Mark Simonson - Proxima Nova Alt Regular-webfont.ttf')
fname='/usr/share/fonts/truetype/Mark Simonson - Proxima Nova Alt Regular-webfont.ttf'
from matplotlib.ticker import FormatStrFormatter
import pylab
import numpy as np
import seaborn as sns
from scipy.integrate import odeint
import itertools
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pandas as pd
from operator import itemgetter
from scipy import stats

#GET KEY FROM DATAFRAME
def get_key(df, key):
    kL=df[key].tolist()
    kM={}
    for k in kL:
        kM[k]=1
    kL=[]
    for k in kM:
        kL.append(k)
    return sorted(kL)

def plot_heatmap(X1,Y1,R,mut,ylabel,vmax,X,Y):
    text=mut.split()
    prot=text[0][:-2]
    fondo=text[1]
    if len(text)>2:
        fondo=fondo+' '+text[2]
    fig=plt.figure()
    ax=plt.gca()
    I=ax.pcolormesh(X,Y,R, vmax=vmax, vmin=0, rasterized=True, cmap='coolwarm')#'RdBu_r')
    ax.set_xscale('log')
    # ax.set_xlim([10**(-2.5),10**(2.5)])
    cr=0.2
    if 'growth' in ylabel: 
        cr=0.4
    C=np.arange(0.1,vmax,cr)
    CS=ax.contour(X,Y,R,C, linewidths=0.3, colors='k')#cmap='binary')
    CL=ax.clabel(CS, C, inline=True, fmt='%.1f', fontsize=7, use_clabeltext=True)
    cbar=fig.colorbar(I, ax=ax)
    cbar.ax.set_yticklabels(cbar.ax.get_yticks(), fontname='Helvetica')
    cbar.ax.set_ylabel(ylabel, fontname='Helvetica', fontsize=15)
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_ylabel('Daylength', fontname='Helvetica', fontsize=15)
    ax.set_xlabel(prot+'  activity', fontname='Helvetica', fontsize=15)
    ax.set_yticks(range(0,25,4))
    ax.set_yticklabels(ax.get_yticks(), fontname='Helvetica')
    ax.set_xticklabels(ax.get_xticks(), fontname='Helvetica')
    ax.set_title(fondo)
    # ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    plt.minorticks_off()
    ax.set_ylim([0,24])
    fig.savefig(dir+'figures/effect_of_'+prot+'_sobre_fondo_'+fondo+'_relative'+ylabel+'_heatmap.pdf', bbox_inches='tight')
    fig.savefig(dir+'figures/effect_of_'+prot+'_sobre_fondo_'+fondo+'_relative'+ylabel+'_heatmap.png', bbox_inches='tight')
    plt.close()

dir='/home/pcatalan/MEGA/plants_ode/'
M=['COP1ox WT', 'COP1ox phyb', 'COP1ox elf3-8', 'COP1ox phyb elf3-8',
	'PHYBox WT', 'PHYBox cop1-4', 'PHYBox elf3-8', 'PHYBox cop1-4 elf3-8',
	'ELF3ox WT', 'ELF3ox phyb', 'ELF3ox cop1-4', 'ELF3ox cop1-4 phyb']
for mut in M:
    print(mut)
    DF=pd.read_csv(dir+'results/heatmap_'+mut+'.csv')
    X=get_key(DF,'Prot')    
    Y=get_key(DF,'D')
    X=np.array(X)
    R=np.zeros((len(Y),len(X)))
    R1=np.zeros((len(Y),len(X)))
    R2=np.zeros((len(Y),len(X)))
    for i,x in enumerate(X):
        for j,y in enumerate(Y):
            A=DF[(DF['Prot']==x) & (DF['D']==y)]
            A22=A[A['T']==22]
            A28=A[A['T']==28]
            R[j,i]=(A28['Growth'].mean()-A22['Growth'].mean())/A22['Growth'].mean()
            R1[j,i]=A28['Growth'].mean()
            R2[j,i]=A22['Growth'].mean()
    # np.savetxt(dir+'results/phyB_matrix.txt', R, fmt='%.3f')
    # np.savetxt(dir+'results/phyB_xaxis.txt', X, fmt='%.6f')
    # np.savetxt(dir+'results/phyB_yaxis.txt', Y, fmt='%.1f')
    yshift=0.5#0.25
    yadd=24+yshift
    Y1=np.array(Y)-yshift
    Y1=list(Y)+[yadd]
    xshift=0.1#0.05
    xadd=10**(3+xshift)
    X1=np.log10(X)-xshift
    X1=list(10**X1)+[xadd]
    # X=np.log10(X)
    
    vmax=R.max()
    print(vmax)
    # continue
    vmax=14
    Vmax=3.0
    plot_heatmap(X1,Y1,R,mut,'Hypocotyl Thermoelongation',Vmax,X,Y)
    # plot_heatmap(X1,Y1,R1,mut,'growth at 28',vmax,X,Y)
    # plot_heatmap(X1,Y1,R2,mut,'growth at 22',vmax,X,Y)
    
    # #CURVES
    # D=[0,4,8,12,16,20,24]
    # vals=np.linspace(0,1,len(D))#each class will have one color
    # cmap=plt.cm.colors.ListedColormap(plt.cm.viridis(vals))
    # fig=plt.figure()
    # ax=plt.gca()
    # text=mut.split()
    # prot=text[0][:-2]
    # fondo=text[1]
    # if len(text)>2:
    #     fondo=fondo+' '+text[2]
    # for i,d in enumerate(D):
    #     A=DF[DF['D']==d]
    #     A22=A[A['T']==22]
    #     A28=A[A['T']==28]
    #     x=np.array(A22['Prot'])
    #     y=np.array(A28['Growth'])-np.array(A22['Growth'])
    #     DFnew=pd.DataFrame({'Prot':x, 'Thermo':y})
    #     DFnew.to_excel(dir+'results/effect_of_'+prot+'_sobre_fondo_'+fondo+'_'+str(d)+'.xls', float_format='%.6f', index=False)
    #     ax.plot(x,y,color=cmap(i/len(D)), label='D='+str(d))
    # ax.legend()
    # ax.set_ylabel('Thermomorphogenetic response', fontproperties=prop, fontsize=15)
    # ax.set_xlabel('phyB activity', fontproperties=prop, fontsize=15)
    # ax.set_xscale('log')
    # fig.savefig(dir+'figures/effect_of_'+prot+'_sobre_fondo_'+fondo+'_curves.pdf', bbox_inches='tight')



# plot_heatmap(X1,Y1,R1,'phyB','growth at 28',vmax,X,Y)
# plot_heatmap(X1,Y1,R2,'phyB','growth at 22',vmax,X,Y)
# #%%
# dir='/home/pcatalan/MEGA/plants_ode/'
# DF=pd.read_csv(dir+'results/effect_of_cop1_01_paper.csv')
# X=get_key(DF,'COP1')
# Y=get_key(DF,'D')
# X=np.array(X)
# R=np.zeros((len(Y),len(X)))
# R1=np.zeros((len(Y),len(X)))
# R2=np.zeros((len(Y),len(X)))
# for i,x in enumerate(X):
#     for j,y in enumerate(Y):
#         A=DF[(DF['COP1']==x) & (DF['D']==y)]
#         A22=A[A['T']==22]
#         A28=A[A['T']==28]
#         R[j,i]=A28['Growth'].mean()-A22['Growth'].mean()
#         R1[j,i]=A28['Growth'].mean()
#         R2[j,i]=A22['Growth'].mean()
# np.savetxt(dir+'results/cop1_matrix.txt', R, fmt='%.3f')
# np.savetxt(dir+'results/cop1_xaxis.txt', X, fmt='%.6f')
# np.savetxt(dir+'results/cop1_yaxis.txt', Y, fmt='%.1f')
# yshift=0.25
# yadd=24+yshift
# Y1=np.array(Y)-yshift
# Y1=list(Y)+[yadd]
# xshift=0.05
# xadd=10**(3+xshift)
# X1=np.log10(X)-xshift
# X1=list(10**X1)+[xadd]
# # X=np.log10(X)

# vmax=max(R1.max(),R2.max())
# plot_heatmap(X1,Y1,R,'COP1','thermomorphogenic response', Vmax,X,Y)
# # plot_heatmap(X1,Y1,R1,'COP1','growth at 28',vmax,X,Y)
# # plot_heatmap(X1,Y1,R2,'COP1','growth at 22',vmax,X,Y)


# dir='/home/pcatalan/MEGA/plants_ode/'
# DF=pd.read_csv(dir+'results/effect_of_elf3_01_paper.csv')
# X=get_key(DF,'ELF3')
# Y=get_key(DF,'D')
# X=np.array(X)
# R=np.zeros((len(Y),len(X)))
# R1=np.zeros((len(Y),len(X)))
# R2=np.zeros((len(Y),len(X)))
# for i,x in enumerate(X):
#     for j,y in enumerate(Y):
#         A=DF[(DF['ELF3']==x) & (DF['D']==y)]
#         A22=A[A['T']==22]
#         A28=A[A['T']==28]
#         R[j,i]=A28['Growth'].mean()-A22['Growth'].mean()
#         R1[j,i]=A28['Growth'].mean()
#         R2[j,i]=A22['Growth'].mean()
# np.savetxt(dir+'results/elf3_matrix.txt', R, fmt='%.3f')
# np.savetxt(dir+'results/elf3_xaxis.txt', X, fmt='%.6f')
# np.savetxt(dir+'results/elf3_yaxis.txt', Y, fmt='%.1f')
# yshift=0.25
# yadd=24+yshift
# Y1=np.array(Y)-yshift
# Y1=list(Y)+[yadd]
# xshift=0.05
# xadd=10**(3+xshift)
# X1=np.log10(X)-xshift
# X1=list(10**X1)+[xadd]
# # X=np.log10(X)
           
# vmax=max(R1.max(),R2.max())
# plot_heatmap(X1,Y1,R,'ELF3','thermomorphogenic response',Vmax,X,Y)
# # plot_heatmap(X1,Y1,R1,'ELF3','growth at 28',vmax,X,Y)
# # plot_heatmap(X1,Y1,R2,'ELF3','growth at 22',vmax,X,Y)


# # np.savetxt(dir+'results/phyb_matrix.txt', R, fmt='%.3f')
# # np.savetxt(dir+'results/phyb_xaxis.txt', X, fmt='%.6f')
# # np.savetxt(dir+'results/phyb_yaxis.txt', Y, fmt='%.1f')


# #%%
# import matplotlib
# import matplotlib.pyplot as plt
# import matplotlib.font_manager as fm
# prop = fm.FontProperties(fname='/usr/share/fonts/truetype/Mark Simonson - Proxima Nova Alt Regular-webfont.ttf')
# fname='/usr/share/fonts/truetype/Mark Simonson - Proxima Nova Alt Regular-webfont.ttf'
# from matplotlib.ticker import FormatStrFormatter
# import pylab
# import numpy as np
# import seaborn as sns
# from scipy.integrate import odeint
# import itertools
# from scipy.optimize import curve_fit
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
# import pandas as pd
# from operator import itemgetter
# from scipy import stats

# #GET KEY FROM DATAFRAME
# def get_key(df, key):
#     kL=df[key].tolist()
#     kM={}
#     for k in kL:
#         kM[k]=1
#     kL=[]
#     for k in kM:
#         kL.append(k)
#     return sorted(kL)


# dir='/home/pcatalan/MEGA/plants_ode/'
# DF=pd.read_csv(dir+'results/effect_of_phyb_01_paper.csv')
# X=get_key(DF,'PHYB')
# Y=get_key(DF,'D')
# X=np.array(X)
# D=[0,4,8,12,16,20,24]
# vals=np.linspace(0,1,len(D))#each class will have one color
# cmap=plt.cm.colors.ListedColormap(plt.cm.viridis(vals))
# fig=plt.figure()
# ax=plt.gca()
# for i,d in enumerate(D):
#     A=DF[DF['D']==d]
#     A22=A[A['T']==22]
#     A28=A[A['T']==28]
#     x=np.array(A22['PHYB'])
#     y=np.array(A28['Growth'])-np.array(A22['Growth'])
#     DFnew=pd.DataFrame({'phyb':x, 'Thermo':y})
#     DFnew.to_excel(dir+'results/effect_of_phyb_'+str(d)+'.xls', float_format='%.6f', index=False)
#     ax.plot(x,y,color=cmap(i/len(D)), label='D='+str(d))
# ax.legend()
# ax.set_ylabel('Thermomorphogenetic response', fontproperties=prop, fontsize=15)
# ax.set_xlabel('phyB activity', fontproperties=prop, fontsize=15)
# ax.set_xscale('log')
# fig.savefig(dir+'figures/effect_of_phyb_curves.pdf', bbox_inches='tight')
# # ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# # ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))