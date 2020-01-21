import numpy as np
import pylab
import subprocess
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from matplotlib import rc
import matplotlib as mpl
import sys

mpl.rcParams['text.usetex'] = 1.0
mpl.rcParams['text.latex.unicode'] = 1.0
mpl.rcParams['xtick.major.pad'] = 10.0
mpl.rcParams['ytick.major.pad'] = 5.0

#############################
## DIHEDRAL ANGLES PHI/PSI ##
#############################
subprocess.call(["rm","-rf","Results/PlotDA"]) 
subprocess.call(["mkdir","Results/PlotDA"]) 

DA = np.genfromtxt('Results/datToPlotDA/DA_all_filenames.dat',dtype='U50',delimiter="\n");

for i in range (len(DA)):
     DAi = np.genfromtxt('Results/datToPlotDA/%s.dat'%(DA[i]))
     ##RAMACHANDRAN PLOTS##
     plt.figure(figsize=(10,10))
     plt.scatter(DAi[:,1],DAi[:,2],s=10,color='w',edgecolor='b')   
     plt.xlabel(r'$\phi\ [\mathrm{degrees} ^{\circ}]$',size=25)
     plt.ylabel(r'$\psi\ [\mathrm{degrees} ^{\circ}]$',size=25)
     plt.title('Ramachandran Plot %s'%(DA[i]),size=20)
     plt.grid()
     plt.xlim(-180,180)
     plt.ylim(-180,180)
     ax = plt.axes()
     ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
     ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
     ax.yaxis.set_major_locator(ticker.MultipleLocator(30))
     ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))
     ax.xaxis.grid(True,'major',linewidth=2)
     ax.xaxis.grid(True,'minor')
     ax.yaxis.grid(True,'major',linewidth=2)
     ax.yaxis.grid(True,'minor')
     ax.axhline(0,linewidth=1.5 ,linestyle='--',color='k')
     ax.axvline(0,linewidth=1.5 ,linestyle='--',color='k')
     ax.tick_params(axis='both',which='major',labelsize=20)
     plt.savefig('Results/PlotDA/DA %s.png'%(DA[i]),format='png',bbox_inches='tight')
     #plt.savefig('Results/PlotDA/DA %s.eps'%(DA[i]),format='eps',bbox_inches='tight')
     #plt.savefig('Results/PlotDA/DA %s.tiff'%(DA[i]),format='tiff',bbox_inches='tight')
     plt.close()

     ##ALL DIHEDRAL ANGLES VS TIME##
     plt.figure(figsize=(25,8))
     psi, = plt.plot(DAi[:,0],DAi[:,1],color='#069af3',markeredgecolor='#069af3',linestyle='none',linewidth=0.3,marker='o',markersize=6.5,label='psi',alpha=0.7)
     phi, = plt.plot(DAi[:,0],DAi[:,2],color='#fb5581',markeredgecolor='r',linestyle='none',linewidth=0.3,marker='s',markersize=6,label='phi',alpha=0.7)
     l=plt.legend([psi,phi],[r'$\phi\mid$',r'$\psi$'],bbox_to_anchor=(0.60,1.16),numpoints=1,ncol=2,markerscale=3,prop={'size':40},borderpad=0.05,handletextpad=0.001,columnspacing=0.1)
     plt.ylabel(r'$\phi,\psi\ [\mathrm{degrees} ^{\circ}]$',size=40.0)
     try:
       chi1, = plt.plot(DAi[:,0],DAi[:,3],color='#21fc0d',markeredgecolor='#21fc0d',linestyle='none',markeredgewidth=1.7,marker='*',markersize=7,label='chi1',alpha=0.7)
       l=plt.legend([psi,phi,chi1],[r'$\phi\mid$',r'$\psi\mid$',r'$\chi_{1}$'],bbox_to_anchor=(0.64,1.16),numpoints=1,ncol=3,markerscale=3,prop={'size':40},borderpad=0.05,handletextpad=0.001,columnspacing=0.1)
       plt.ylabel(r'$\phi,\psi,\chi_{1}\ [\mathrm{degrees} ^{\circ}]$',size=40.0)
     except IndexError:
       print ('Some index does not exist')
     try:
       chi2, = plt.plot(DAi[:,0],DAi[:,4],color='#bc13fe',linestyle='none',markeredgewidth=1.7,marker='+',markersize=8,label='chi2',alpha=0.7)
       l=plt.legend([psi,phi,chi1,chi2],[r'$\phi\mid$',r'$\psi\mid$',r'$\chi_{1}\mid$',r'$\chi_{2}$'],bbox_to_anchor=(0.69,1.16),numpoints=1,ncol=4,markerscale=3,prop={'size':40},borderpad=0.05,handletextpad=0.001,columnspacing=0.1)
       plt.ylabel(r'$\phi,\psi,\chi_{1},\chi_{2}\ [\mathrm{degrees} ^{\circ}]$',size=40.0)
     except IndexError:
       print ('Some index does not exist')
     try:
       chi3, = plt.plot(DAi[:,0],DAi[:,5],color='#4b5d16',linestyle='none',markeredgewidth=1.7,marker='x',markersize=7,label='chi3',alpha=0.7)
       l=plt.legend([psi,phi,chi1,chi2,chi3],[r'$\phi\mid$',r'$\psi\mid$',r'$\chi_{1}\mid$',r'$\chi_{2}\mid$',r'$\chi_{3}$'],bbox_to_anchor=(0.74,1.16),numpoints=1,ncol=5,markerscale=3,prop={'size':40},borderpad=0.05,handletextpad=0.001,columnspacing=0.1)
       plt.ylabel(r'$\phi,\psi,\chi_{1},\chi_{2},\chi_{3}\ [\mathrm{degrees} ^{\circ}]$',size=40.0)
     except IndexError:
       print ('Some index does not exist')
     try:
       chi4, = plt.plot(DAi[:,0],DAi[:,6],color='darkgoldenrod',markeredgecolor='darkgoldenrod',linestyle='none',markeredgewidth=3,marker='_',markersize=7.5,label='chi4',alpha=0.7)
       l=plt.legend([psi,phi,chi1,chi2,chi3,chi4],[r'$\phi\mid$',r'$\psi\mid$',r'$\chi_{1}\mid$',r'$\chi_{2}\mid$',r'$\chi_{3}\mid$',r'$\chi_{4}\mid$'],bbox_to_anchor=(0.77,1.16),numpoints=1,ncol=6,markerscale=3,prop={'size':40},borderpad=0.05,handletextpad=0.001,columnspacing=0.1)
       plt.ylabel(r'$\phi,\psi,\chi_{1},\chi_{2},\chi_{3},\chi_{4}\ [\mathrm{degrees} ^{\circ}]$',size=40.0)
     except IndexError:
       print ('Some index does not exist')
     try:
       chi5, = plt.plot(DAi[:,0],DAi[:,7],color='c',markeredgecolor='c',linestyle='none',markeredgewidth=3,marker='|',markersize=7.5,label='chi5',alpha=0.7)
       l=plt.legend([psi,phi,chi1,chi2,chi3,chi4,chi5],[r'$\phi\mid$',r'$\psi\mid$',r'$\chi_{1}\mid$',r'$\chi_{2}\mid$',r'$\chi_{3}\mid$',r'$\chi_{4}\mid$',r'$\chi_{5}$'],bbox_to_anchor=(0.82,1.16),numpoints=1,ncol=7,markerscale=3,prop={'size':40},borderpad=0.05,handletextpad=0.001,columnspacing=0.1)
       plt.ylabel(r'$\phi,\psi,\chi_{1},\chi_{2},\chi_{3},\chi_{4},\chi_{5}\ [\mathrm{degrees} ^{\circ}]$',size=40.0)
     except IndexError:
       print ('Some index does not exist')

     plt.xlim([-10,np.max(DAi[:,0])])
     plt.ylim(-190,190)

     plt.xlabel('Time 'r'$[\mathrm{ns}]$',size=40.0,fontsize=40)
     plt.title('%s'%(DA[i]),x=0.92,y=1.04,size=(40.0))
     ax = plt.axes()
     ax.yaxis.set_major_locator(ticker.MultipleLocator(30))
     ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))
     ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
     ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
     ax.tick_params(axis='both',which='major',labelsize=45)
     #for label in ax.yaxis.get_ticklabels()[1::2]:
      # label.set_visible(False)
     #for label in ax.xaxis.get_ticklabels()[0::2]:
      # label.set_visible(False)

     ax.xaxis.grid(True,'major',linewidth=2)
     ax.xaxis.grid(True,'minor',linewidth=1)
     ax.yaxis.grid(True,'major',linewidth=2)
     ax.yaxis.grid(True,'minor',linewidth=1)
     ax.axhline(0,linewidth=1 ,linestyle='--',color='k')
     plt.savefig('Results/PlotDA/DA time plots %s.png'%(DA[i]),format='png',bbox_inches='tight')
     #plt.savefig('Results/PlotDA/DA time plots %s.tiff'%(DA[i]),format='tiff',bbox_inches='tight')
     plt.close("all")

print ('Finished DA plots')
