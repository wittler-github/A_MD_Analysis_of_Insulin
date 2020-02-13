import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pylab
import subprocess
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from matplotlib import rc
import matplotlib as mpl
import sys
import matplotlib.gridspec as gridspec

#mpl.rc('font', family='serif', serif='cm10')
mpl.rcParams['text.usetex'] = 1.0
#mpl.rcParams['text.latex.preamble'] = [r'\boldmath']
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
mpl.rcParams['axes.axisbelow'] = True
#mpl.rcParams["font"] = "sansserif"
#mpl.rcParams["axes.labelweight"] = "bold"

subprocess.call(["rm","-rf","Results/PlotCSNP"]) 
subprocess.call(["mkdir","Results/PlotCSNP"]) 

CSNP_filenames = np.genfromtxt('Results/datToPlotCSNP/CSNP_all_filenames.dat',dtype='U50',delimiter="\t",autostrip=True);
QFCstat = np.genfromtxt('Results/datToPlotCSNP/QFCstat_all_filenames.dat',dtype=('U50',None,None,None,None,None,None,None,None),delimiter="\t");
SASA = np.genfromtxt('Results/datToPlotCSNP/SASA(MD).dat')

###################################################################
## Calculables  ##
###################################################################
for i in range(5):#len(CSNP_filenames)):
  if (CSNP_filenames[i] != QFCstat[i][0]):
    print (' NOT EQUAL '+'CSNP_filenames['+'%i'%(i)+']'+' != '+'QFCstat['+'%i'%(i)+'][0] :'+' %s'%(CSNP_filenames[i])+' != '+' %s'%(QFCstat[i][0])+' "Exiting Script"')
    sys.exit()
  CSNPi = np.genfromtxt('Results/datToPlotCSNP/%s.dat'%(CSNP_filenames[i]))
  fig=plt.figure(figsize=(30.0,10.0))
  plt.xlim([0,np.max(CSNPi[:,0])])
  plt.rcParams['xtick.labelsize'] = 55.0
  plt.rcParams['ytick.labelsize'] = 55.0

  ax = plt.gca()
  ax2 = ax.twinx()
  axl, = ax.plot(SASA[:,0],SASA[:,1],color='darkolivegreen',markeredgecolor='darkolivegreen',markeredgewidth=2.5,linewidth=3.0,linestyle='-',marker='',markersize=7.0,alpha=1.0,zorder=6)
  ax2l, = ax2.plot(CSNPi[:,0],CSNPi[:,1],color='crimson',markeredgecolor='crimson',markeredgewidth=1.7,linewidth=3.0,linestyle='-',marker='',markersize=7.0,alpha=0.5,zorder=2)
  l1 = ax.legend([axl],[r'$\mid$'+r'$\langle x \rangle$=3986'+r' $\mathrm{\AA}^2\mid$'+r' $SD_{x_{i}}$=179.1'+r' $\mathrm{\AA}^2\mid$'],bbox_to_anchor=(0.51,1.15),numpoints=1.0,ncol=1.0,markerscale=4.0,prop={'size':50.0},borderpad=0.2,framealpha=0.4)
  l2 = ax2.legend([ax2l],[r'$\mid$'+r'$\langle x \rangle$='+'%.1f'%(QFCstat[i][1])+r' $\mathrm{\AA}\mid$'+r' $SD_{x_{i}}$='+'%.1f'%(QFCstat[i][2])+r' $\mathrm{\AA}\mid$'],bbox_to_anchor=(1.01,1.15),numpoints=1.0,ncol=1.0,markerscale=4.0,prop={'size':50.0},borderpad=0.35,framealpha=0.4)

  ax.set_xlim([0.0,np.max(SASA[:,0])])
  ax.set_ylim([np.min(SASA[:,1])-(np.min(SASA[:,1])/30),np.max(SASA[:,1])+(np.max(SASA[:,1])/15)])
  ax2.set_ylim([np.min(CSNPi[:,1])-(np.min(CSNPi[:,1])/30),np.max(CSNPi[:,1])+(np.max(CSNPi[:,1])/15)])

  plt.xlabel('Time [ns]',size=50.0)
  ax.set_ylabel('SASA [Å$^2$]',size=55.0,color='darkolivegreen')
  ax2.set_ylabel('Nr Water(O) within '+'%.0f'%(QFCstat[i][7])+' Å',size=55.0,color='crimson')
  ax.yaxis.set_label_position("left")
  ax2.yaxis.set_label_position("right")
  plt.rc('font',weight='bold')
  ax.xaxis.set_major_locator(ticker.MultipleLocator(100.0))
  ax.xaxis.set_minor_locator(ticker.MultipleLocator(10.0))
  #ax.yaxis.set_major_locator(ticker.MultipleLocator(100.0))
  #ax.yaxis.set_minor_locator(ticker.MultipleLocator(10.0))
  #ax2.yaxis.set_major_locator(ticker.MultipleLocator(100.0))
  #ax2.yaxis.set_minor_locator(ticker.MultipleLocator(10.0))
  ax.set_axisbelow(True)
  ax.xaxis.grid(True,'major',color='dimgray',linewidth=3.0)
  ax.xaxis.grid(True,'minor',color='dimgray',linewidth=1.5)
  ax2.set_axisbelow(True)
  #ax.yaxis.grid(True,'major',color='green',linewidth=3.0)
  #ax2.yaxis.grid(True,'major',color='red',linewidth=3.0, alpha=0.5)
  #ax.yaxis.grid(False,'minor',color='green',linewidth=1.5,zorder=4)
  #ax2.yaxis.grid(False,'minor',color='red',linewidth=1.5,zorder=5)
  #ax1.plot(CSNPi[:,0],CSNPi[:,1],color='crimson',linestyle='-',linewidth=0.8)
  #ax2.plot(CSNPi[:,0],CSNPi[:,2],color='blue',linestyle='-',linewidth=0.8)
  #ax2.plot(CSNPi[:,0],CSNPi[:,3],color='blue',linestyle='-',linewidth=0.8)

  plt.savefig('Results/PlotCSNP/%s.png'%(CSNP_filenames[i]),format='png',bbox_inches='tight')
  #plt.savefig('Results/PlotCSNP/%s.tiff'%(CSNP_filenames[i][0]),format='tiff',bbox_inches='tight')
  plt.close()
