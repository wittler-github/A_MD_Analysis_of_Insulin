import numpy as np
import pylab
import subprocess
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from matplotlib import rc
import matplotlib as mpl
import sys
import matplotlib.gridspec as gridspec

mpl.rcParams['text.usetex'] = 1.0
#mpl.rcParams['text.latex.unicode'] = 1.0
#mpl.rcParams['xtick.major.pad'] = 10.0
#mpl.rcParams['ytick.major.pad'] = 5.0

mpl.rcParams['font.family'] = "serif"
mpl.rcParams['font.serif'] = "Times New Roman"
mpl.rcParams['text.latex.unicode'] = 1.0
mpl.rcParams['xtick.major.pad'] = 10.0
mpl.rcParams['ytick.major.pad'] = 5.0

mpl.rcParams['xtick.major.size'] = 6.0
mpl.rcParams['xtick.major.width'] = 3.0
mpl.rcParams['xtick.minor.size'] = 4.0
mpl.rcParams['xtick.minor.width'] = 2.0
mpl.rcParams['ytick.major.size'] = 6.0
mpl.rcParams['ytick.major.width'] = 3.0
mpl.rcParams['ytick.minor.size'] = 4.0
mpl.rcParams['ytick.minor.width'] = 2.0
mpl.rcParams['xtick.color'] = "black"
mpl.rcParams['ytick.color'] = "black"
mpl.rcParams['axes.labelcolor'] = "black"
mpl.rcParams['grid.alpha'] = 0.7

other_size=32
label_size=33
#############################
## DIHEDRAL ANGLES PHI/PSI ##
#############################
subprocess.call(["rm","-rf","Results/PlotDA"]) 
subprocess.call(["mkdir","Results/PlotDA"]) 

#DAi = np.genfromtxt('Results/datToPlotDA/AllresidDA.dat',delimiter=' ',usecols=np.arange(0,5));
DAi = np.genfromtxt('../Results/datToPlotDA/AllresidDA.dat');
E_DAi = np.genfromtxt('Results/datToPlotDA/E_AllresidDA.dat');
print(E_DAi)
E_DAi_stat = np.genfromtxt('Results/datToPlotDA/E_QFC_DA_stat_all_filenames.dat');
#,dtype=(int,None,None,None,None,None,None,None),delimiter='\n',usemask=True

t_DAi=DAi
t_E_DAi=E_DAi
t_array=E_DAi
t_array[:,1:7]=np.where((t_E_DAi[:,1:7]==1.0)&(t_E_DAi[:,1:7]!=0.0),t_DAi[:,1:7],np.nan)
Ep_DAi=t_array
print(E_DAi)
##ALL DIHEDRAL ANGLES VS TIME##
plt.figure(figsize=(25,11))
gs = gridspec.GridSpec(2, 1,height_ratios=[1.5,0.5])
ax=plt.subplot(gs[0])
E_ax=plt.subplot(gs[1])

#ax = plt.gca()
#E_ax = ax.twinx()

psi, = ax.plot(DAi[:,0],DAi[:,1],color='#069af3',markeredgecolor='#069af3',linestyle='none',linewidth=0.3,marker='o',markersize=11,label='psi',alpha=0.7)
phi, = ax.plot(DAi[:,0],DAi[:,2],color='#fb5581',markeredgecolor='r',linestyle='none',linewidth=0.3,marker='s',markersize=10,label='phi',alpha=0.7)
l=ax.legend([psi,phi],[r'$\phi\mid$',r'$\psi$'],bbox_to_anchor=(0.60,1.10),numpoints=1,ncol=2,markerscale=1.5,prop={'size':other_size},borderpad=0.05,handletextpad=0.001,columnspacing=0.1)
ax.set_ylabel(r'$\phi,\psi\ [\mathrm{degrees} ^{\circ}]$',size=other_size)

Ep_phi, = ax.plot(Ep_DAi[:,0],Ep_DAi[:,1],color='w',markeredgecolor='k',linestyle='none',markeredgewidth=2.0,marker='H',markersize=15,label='Ep_phi',alpha=0.7,zorder=0)
Ep_psi, = ax.plot(Ep_DAi[:,0],Ep_DAi[:,2],color='w',markeredgecolor='k',linestyle='none',markeredgewidth=2.0,marker='H',markersize=15,label='Ep_psi',alpha=0.7,zorder=0)

Es_phi, = E_ax.plot(E_DAi_stat[:,0],E_DAi_stat[:,1],color='#069af3',markeredgecolor='#069af3',linestyle='none',linewidth=0.3,marker='o',markersize=10,label='Es_phi',alpha=0.8)
Es_psi, = E_ax.plot(E_DAi_stat[:,0],E_DAi_stat[:,2],color='#fb5581',markeredgecolor='#fb5581',linestyle='none',linewidth=0.3,marker='s',markersize=11,label='Es_psi',alpha=0.8)
#Es_l=ax.legend([Es_psi],[r''],numpoints=1,ncol=2,markerscale=1.5,prop={'size':40})
E_ax.set_ylabel('Fraction of time',size=other_size)
try:
  chi1, = ax.plot(DAi[:,0],DAi[:,3],color='#21fc0d',markeredgecolor='#21fc0d',linestyle='none',markeredgewidth=1.7,marker='*',markersize=12,label='chi1',alpha=0.7)
  l=ax.legend([psi,phi,chi1],[r'$\phi\mid$',r'$\psi\mid$',r'$\chi_{1}$'],bbox_to_anchor=(0.64,1.10),numpoints=1,ncol=3,markerscale=1.5,prop={'size':other_size},borderpad=0.05,handletextpad=0.001,columnspacing=0.1)
  ax.set_ylabel(r'$\phi,\psi,\chi_{1}\ [\mathrm{degrees} ^{\circ}]$',size=other_size)
  Ep_chi1, = ax.plot(Ep_DAi[:,0],Ep_DAi[:,3],color='w',markeredgecolor='k',linestyle='none',markeredgewidth=2.0,marker='H',markersize=15,label='Ep_chi1',alpha=0.7,zorder=0)
  Es_chi1, = E_ax.plot(E_DAi_stat[:,0],E_DAi_stat[:,3],color='#21fc0d',markeredgecolor='#21fc0d',linestyle='none',markeredgewidth=1.7,marker='*',markersize=12,label='Es_chi1',alpha=0.8)
except IndexError:
  print ('Some index does not exist')
try:
  chi2, = ax.plot(DAi[:,0],DAi[:,4],color='#bc13fe',linestyle='none',markeredgewidth=3.5,marker='+',markersize=11,label='chi2',alpha=0.7)
  l=ax.legend([psi,phi,chi1,chi2],[r'$\phi\mid$',r'$\psi\mid$',r'$\chi_{1}\mid$',r'$\chi_{2}$'],bbox_to_anchor=(0.69,1.10),numpoints=1,ncol=4,markerscale=1.5,prop={'size':other_size},borderpad=0.05,handletextpad=0.001,columnspacing=0.1)
  ax.set_ylabel(r'$\phi,\psi,\chi_{1},\chi_{2}\ [\mathrm{degrees} ^{\circ}]$',size=other_size)
  Ep_chi2, = ax.plot(Ep_DAi[:,0],Ep_DAi[:,4],color='w',markeredgecolor='k',linestyle='none',markeredgewidth=2.0,marker='H',markersize=15,label='Ep_chi2',alpha=0.7,zorder=0)
  Es_chi2, = E_ax.plot(E_DAi_stat[:,0],E_DAi_stat[:,4],color='#bc13fe',linestyle='none',markeredgewidth=3.5,marker='+',markersize=11,label='Es_chi2',alpha=0.8)
except IndexError:
  print ('Some index does not exist')
try:
  chi3, = ax.plot(DAi[:,0],DAi[:,5],color='#4b5d16',linestyle='none',markeredgewidth=3.5,marker='x',markersize=9,label='chi3',alpha=0.7)
  l=ax.legend([psi,phi,chi1,chi2,chi3],[r'$\phi\mid$',r'$\psi\mid$',r'$\chi_{1}\mid$',r'$\chi_{2}\mid$',r'$\chi_{3}$'],bbox_to_anchor=(0.74,1.10),numpoints=1,ncol=5,markerscale=1.5,prop={'size':40},borderpad=0.05,handletextpad=0.001,columnspacing=0.1)
  ax.set_ylabel(r'$\phi,\psi,\chi_{1},\chi_{2},\chi_{3}\ [\mathrm{degrees} ^{\circ}]$',size=other_size)
  Ep_chi3, = ax.plot(Ep_DAi[:,0],Ep_DAi[:,5],color='w',markeredgecolor='k',linestyle='none',markeredgewidth=2.0,marker='H',markersize=15,label='Ep_chi3',alpha=0.7,zorder=0)
  Es_chi3, = E_ax.plot(E_DAi_stat[:,0],E_DAi_stat[:,5],color='#4b5d16',linestyle='none',markeredgewidth=3.5,marker='x',markersize=9,label='Es_chi3',alpha=0.8)
except IndexError:
  print ('Some index does not exist')
try:
  chi4, = ax.plot(DAi[:,0],DAi[:,6],color='darkgoldenrod',markeredgecolor='darkgoldenrod',linestyle='none',markeredgewidth=3.5,marker='_',markersize=12,label='chi4',alpha=0.7)
  l=ax.legend([psi,phi,chi1,chi2,chi3,chi4],[r'$\phi\mid$',r'$\psi\mid$',r'$\chi_{1}\mid$',r'$\chi_{2}\mid$',r'$\chi_{3}\mid$',r'$\chi_{4}\mid$'],bbox_to_anchor=(0.77,1.10),numpoints=1,ncol=6,markerscale=1.5,prop={'size':other_size},borderpad=0.05,handletextpad=0.001,columnspacing=0.1)
  ax.set_ylabel(r'$\phi,\psi,\chi_{1},\chi_{2},\chi_{3},\chi_{4}\ [\mathrm{degrees} ^{\circ}]$',size=other_size)
  Ep_chi4, = ax.plot(Ep_DAi[:,0],Ep_DAi[:,6],color='w',markeredgecolor='k',linestyle='none',markeredgewidth=2.0,marker='H',markersize=15,label='Ep_psi',alpha=0.7,zorder=0)
  Es_chi4, = E_ax.plot(E_DAi_stat[:,0],E_DAi_stat[:,6],color='darkgoldenrod',markeredgecolor='darkgoldenrod',linestyle='none',markeredgewidth=3.5,marker='_',markersize=12,label='Es_chi4',alpha=0.8)
except IndexError:
  print ('Some index does not exist')
try:
  chi5, = ax.plot(DAi[:,0],DAi[:,7],color='c',markeredgecolor='c',linestyle='none',markeredgewidth=3.5,marker='|',markersize=12,label='chi5',alpha=0.7)
  l=ax.legend([psi,phi,chi1,chi2,chi3,chi4,chi5],[r'$\phi\mid$',r'$\psi\mid$',r'$\chi_{1}\mid$',r'$\chi_{2}\mid$',r'$\chi_{3}\mid$',r'$\chi_{4}\mid$',r'$\chi_{5}$'],bbox_to_anchor=(0.82,1.10),numpoints=1,ncol=7,markerscale=1.5,prop={'size':other_size},borderpad=0.05,handletextpad=0.001,columnspacing=0.1)
  ax.set_ylabel(r'$\phi,\psi,\chi_{1},\chi_{2},\chi_{3},\chi_{4},\chi_{5}\ [\mathrm{degrees} ^{\circ}]$',size=other_size)
  Ep_chi5, = ax.plot(Ep_DAi[:,0],Ep_DAi[:,7],color='w',markeredgecolor='k',linestyle='none',markeredgewidth=2.0,marker='H',markersize=15,label='Ep_chi5',alpha=0.7,zorder=0)
  Es_chi5, = E_ax.plot(E_DAi_stat[:,0],E_DAi_stat[:,7],color='c',markeredgecolor='c',linestyle='none',markeredgewidth=3.5,marker='|',markersize=12,label='Es_chi5',alpha=0.8)
except IndexError:
  print ('Some index does not exist')

ax.set_xlim([0.5,51.5])
ax.set_ylim([-190,190])
E_ax.set_ylim([-0.05,1.05])
E_ax.set_xlim([0.5,51.5])
ax.yaxis.tick_left()
E_ax.set_xticklabels([])
E_ax.yaxis.tick_left()
ax.yaxis.set_label_position("left")
E_ax.yaxis.set_label_position("left")
plt.tight_layout()

E_ax.set_xlabel('Residue nr of A-chain/B-chain',size=other_size,fontsize=40)
#plt.title('%s'%(),x=0.92,y=1.04,size=(40.0))
#axes = ax.axes()
#E_axes = ax.axes()
ax.yaxis.set_major_locator(ticker.MultipleLocator(30))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))
ax.xaxis.set_major_locator(ticker.MultipleLocator(3))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
E_ax.yaxis.set_major_locator(ticker.MultipleLocator(0.25))
#E_ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
E_ax.xaxis.set_major_locator(ticker.MultipleLocator(3))
E_ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))

ax.tick_params(axis='both',which='major',labelsize=label_size)
E_ax.tick_params(axis='both',which='major',labelsize=label_size)
#for label in ax.yaxis.get_ticklabels()[1::2]:
 # label.set_visible(False)
#for label in ax.xaxis.get_ticklabels()[0::2]:
 # label.set_visible(False)

ax.xaxis.grid(True,'major',linewidth=2)
ax.xaxis.grid(True,'minor',linewidth=1)
ax.yaxis.grid(True,'major',linewidth=2)
ax.yaxis.grid(True,'minor',linewidth=1)
E_ax.xaxis.grid(True,'major',linewidth=2)
E_ax.xaxis.grid(True,'minor',linewidth=1)
E_ax.yaxis.grid(True,'major',linewidth=2)
E_ax.yaxis.grid(False,'minor',linewidth=1)
#E_ax.axhline(0,linewidth=1 ,linestyle='--',color='k')
#plt.savefig('Results/PlotDA/DA time plots MS_AllresidDA.png',format='png',bbox_inches='tight')
plt.savefig('Results/PlotDA/DA time plots MS_AllresidDA.tiff',format='tiff',bbox_inches='tight',dpi=200)
plt.savefig('Results/PlotDA/DA time plots MS_AllresidDA.svg',format='svg',bbox_inches='tight')
#plt.savefig('Results/PlotDA/DA time plots %s.tiff'%(DA[i]),format='tiff',bbox_inches='tight')
plt.clf()
plt.cla()
plt.close("all")

print ('Finished DA plots')
