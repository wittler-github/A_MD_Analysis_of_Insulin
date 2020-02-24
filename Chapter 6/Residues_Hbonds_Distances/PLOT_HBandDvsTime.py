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
print ("MPL version ")
print (matplotlib.__version__)

mpl.rcParams['text.usetex'] = 1.0
mpl.rcParams['font.family'] = "serif"
mpl.rcParams['font.serif'] = "Times New Roman"
mpl.rcParams['text.latex.unicode'] = 1.0
mpl.rcParams['xtick.major.pad'] = 10.0
mpl.rcParams['ytick.major.pad'] = 5.0
mpl.rcParams['xtick.major.size'] = 6.0
mpl.rcParams['xtick.major.width'] = 3.0
mpl.rcParams['xtick.minor.size'] = 3.0
mpl.rcParams['xtick.minor.width'] = 1.5
mpl.rcParams['ytick.major.size'] = 6.0
mpl.rcParams['ytick.major.width'] = 3.0
mpl.rcParams['ytick.minor.size'] = 3.0
mpl.rcParams['ytick.minor.width'] = 1.5
mpl.rcParams['xtick.color'] = "black"
mpl.rcParams['ytick.color'] = "black"
mpl.rcParams['axes.labelcolor'] = "black"
mpl.rcParams['grid.alpha'] = 0.7

label_size=30
Llabel_size=27

'''
####################
## HYDROGEN BONDS ##
####################
print ('Commencing HB plots')

subprocess.call(["rm","-rf","Results/PlotHBvsTime"]) 
subprocess.call(["mkdir","Results/PlotHBvsTime"]) 

HB = np.genfromtxt('Results/datToPlotHB/HB_all_filenames.dat',dtype='U50',delimiter="\n",skip_header=0,skip_footer=0);
HBstat = np.genfromtxt('Results/datToPlotHB/QFCstat_all_filenames.dat',dtype=('U50',None,None,None,None),delimiter="\t",usecols=np.arange(0,5),skip_header=0,skip_footer=0);
HBperc = np.genfromtxt('Results/datToPlotHB/HBperc_all_filenames.dat',dtype=('U50',None,None,None,None),delimiter="\t",skip_header=0,skip_footer=0);
lenHB = len(HB)-1
lenHBperc = len(HBperc)
    
for i in list(range(54,56)):
  if (HB[i] != HBstat[i][0]) | (HB[i] != HBperc[i][0]):
    print (' NOT EQUAL '+'HB['+'%i'%(i)+']'+' != '+'HBstat['+'%i'%(i)+'][0] :'+' %s'%(HB[i])+' != '+' %s'%(HBstat[i][0])+' "Exiting Script"')
    sys.exit()
  if (HB[i] == HBstat[i][0]) & (HB[i] == HBperc[i][0]) & (HBperc[i][2]>=5):
    print (' EQUAL '+'HB['+'%i'%(i)+']'+' = '+'HBstat['+'%i'%(i)+'][0] & '+'HBperc['+'%i'%(i)+'][0] =  :'+' %s'%(HB[i])+' = '+' %s'%(HBstat[i][0])+' & %s'%(HBperc[i][0]))
    HBi = np.genfromtxt('Results/datToPlotHB/%s.dat'%(HB[i]))
    plt.figure(figsize=(24.0,8.0))
    #HBiSX = (0.001*HBi[:,0])
    HBiSX = (HBi[:,0])
    ax = plt.gca()
    ax2 = ax.twinx()
    axl, = ax.plot(HBiSX[::5],HBi[:,1][::5],color='m',markeredgecolor='m',markeredgewidth=4.0,linewidth=1.0,linestyle='none',marker='|',markersize=8,alpha=0.7,zorder=2)
    ax2l, = ax2.plot(HBiSX[::5],HBi[:,2][::5],color='g',markeredgecolor='g',markeredgewidth=2,linewidth=1.0,linestyle='none',marker='^',markersize=7.0,alpha=0.7,zorder=2)
    thirtee = a=np.empty(len(HBi[:,0])); a.fill(30) 
    threepfive = a=np.empty(len(HBi[:,0])); a.fill(3.5)
    ax.plot(HBiSX,threepfive,color='m',markeredgecolor='m',linestyle='-',marker='.',zorder=0)
    ax2.plot(HBiSX,thirtee,color='g',markeredgecolor='g',linestyle='-',marker='.',zorder=0)
    ax.set_xlim([0.0,np.max(HBiSX)])
    ax.set_ylim([-0.5,10.5])
    ax2.set_ylim([-10.0,190.0])
    ax.yaxis.tick_left()
    ax2.yaxis.tick_right()
    ax.yaxis.set_label_position("left")
    ax2.yaxis.set_label_position("right")
    ax.tick_params(axis='x',which='major',labelsize=label_size, pad=0.2)
    ax.tick_params(axis='y',which='mahor',labelsize=label_size,right=0,labelleft='on')
    ax2.tick_params(axis='y',which='major',labelsize=label_size,right=1,labelright='on')
    plt.rcParams['xtick.labelsize'] = label_size
    plt.rcParams['ytick.labelsize'] = label_size
    plt.rc('font',weight='bold')

    l1 = ax.legend([axl],[r'$r_{AD}[\mathrm{\AA}]\mid$'+r'$\langle x \rangle$='+'%.1f'%(HBstat[i][1])+r'$\mathrm{\AA}\mid$'+r' $SD_{x_{i}}$='+'%.1f'%(HBstat[i][2])+r'$\mathrm{\AA}\mid$'],bbox_to_anchor=(0.19,1.14),numpoints=1.0,ncol=1.0,markerscale=2.0,prop={'size':Llabel_size},borderpad=0.4,frameon=1,edgecolor='k',framealpha=1.0)
    l2 = ax2.legend([ax2l],[r'$\langle x \rangle$='+'%.1f'%(HBstat[i][3])+r'$^{\circ}\mid$'+r' $SD_{x_{i}}$='+'%.1f'%(HBstat[i][4])+r'$^{\circ}\mid$'+r'        $\varphi [\mathrm{degrees}^{\circ}]$'],bbox_to_anchor=(0.999,1.14),numpoints=1.0,ncol=1.0,markerscale=2.0,prop={'size':Llabel_size},borderpad=0.41,frameon=1,edgecolor='k',framealpha=1.0,markerfirst=False)
    
    ax.text(-33,-1.8,'%s'%(HB[i]),size=(label_size),rotation=0,bbox=dict(facecolor='w',edgecolor='black'))
    ax.text(700,-1.8,'Time 'r'$[\mathrm{ns}]$',size=label_size,fontsize=label_size)

    ax.text(461.0,9.0,r'$( r_{AD} \ < \ 3.5 \dot A ) \ \& \ ( \varphi \ < \ 30^{\circ} ) \ : \ $'+'%.1f'%(HBperc[i][2])+r'$\%$'+r'$\%$',fontsize=41.0,color='k',bbox=dict(facecolor='none',edgecolor='slateblue',boxstyle='Darrow,pad=0.3'))
    ax.text(20.5,9.0,r'$(r_{AD} \ < \ 3.5 \dot A) \ : \ $'+'%.1f'%(HBperc[i][3])+r'$\%$',fontsize=41.0,color='k',bbox=dict(facecolor='none',edgecolor='magenta',boxstyle='sawtooth,pad=0.3'))
    ax.text(1165.0,9.0,r'$(\varphi \ < \ 30^{\circ} ) \ : \ $'+'%.1f'%(HBperc[i][4])+r'$\%$',fontsize=41.0,color='k',bbox=dict(facecolor='none',edgecolor='green',boxstyle='sawtooth,pad=0.3'))

    for tic in ax.get_yticklabels():
      tic.set_color('m')
    for tic in ax2.get_yticklabels():
      tic.set_color('darkgreen')

    ax.xaxis.set_major_locator(ticker.MultipleLocator(100.0))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(10.0))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(30.0))
    ax2.yaxis.set_minor_locator(ticker.MultipleLocator(10.0))

    ax.axhline(0,linestyle='-', color='m',linewidth=2,zorder=0)
    ax.axhline(1,linestyle='-', color='m',linewidth=2,zorder=0)
    ax.axhline(2,linestyle='-', color='m',linewidth=2,zorder=0)
    ax.axhline(3,linestyle='-', color='m',linewidth=2,zorder=0)
    ax.axhline(4,linestyle='-', color='m',linewidth=2,zorder=0)
    ax.axhline(5,linestyle='-', color='m',linewidth=2,zorder=0)
    ax.axhline(6,linestyle='-', color='m',linewidth=2,zorder=0)
    ax.axhline(7,linestyle='-', color='m',linewidth=2,zorder=0)
    ax.axhline(8,linestyle='-', color='m',linewidth=2,zorder=0)
    ax.axhline(9,linestyle='-', color='m',linewidth=2,zorder=0)
    ax.axhline(10,linestyle='-', color='m',linewidth=2,zorder=0)

    ax2.axhline(0,linestyle='-', color='g',linewidth=2,zorder=0)
    ax2.axhline(30,linestyle='-', color='g',linewidth=2,zorder=0)
    ax2.axhline(60,linestyle='-', color='g',linewidth=2,zorder=0)
    ax2.axhline(90,linestyle='-', color='g',linewidth=2,zorder=0)
    ax2.axhline(120,linestyle='-', color='g',linewidth=2,zorder=0)
    ax2.axhline(150,linestyle='-', color='g',linewidth=2,zorder=0)
    ax2.axhline(180,linestyle='-', color='g',linewidth=2,zorder=0)

    ax.xaxis.grid(True,'major',linewidth=2.0, color='k')
    ax.xaxis.grid(True,'minor',color='k')
    ax.set_position([0.1,0.1,0.5,0.8])
    plt.tight_layout()
    plt.savefig('Results/PlotHBvsTime/HBvsTime%s.svg'%(HB[i].replace(" ","")),format='svg',bbox_inches='tight',dpi=None,facecolor='#a9f971')
    #plt.savefig('Results/PlotHBvsTime/HBvsTime%s.jpg'%(HB[i].replace(" ","")),format='jpg',bbox_inches='tight',dpi=None,facecolor='#a9f971')
    plt.close()

print ('Finished HB plots')
'''

#PLOT CA AND GC ON SAME PLOT!!!!!
#########################################################
## CARBON ALPHA DISTANCES & GEOMETRIC CENTER DISTANCES ##
#########################################################
print ('Commencing CA/GC plots')


subprocess.call(["rm","-rf","Results/PlotCAGCvsTime"]) 
subprocess.call(["mkdir","Results/PlotCAGCvsTime"]) 

CAGC = np.genfromtxt('Results/datToPlotCAGC/CAGC_all_filenames.dat',dtype='U50',delimiter="\n");
CAGCstat = np.genfromtxt('Results/datToPlotCAGC/QFCstat_all_filenames.dat',dtype=('U50',None,None,None,None,None,None,None,None,None,None,None,None),delimiter="\t",usecols=np.arange(0,12));

for i in range(2):
  if (CAGC[i] != CAGCstat[i][0]):
    print (' NOT EQUAL '+'(CAGC['+'%i'%(i)+']'+' != '+'CAGCstat['+'%i'%(i)+'][0] :'+' %s'%(CAGC[i])+' != '+' %s)'%(CAGCstat[i][0])+' "Exiting Script"')
    sys.exit()
  if (CAGC[i] == CAGCstat[i][0]):
    print (' EQUAL '+'(CAGC['+'%i'%(i)+']'+' = '+'CAGCstat['+'%i'%(i)+'][0] :'+' %s'%(CAGC[i])+' = '+' %s)'%(CAGCstat[i][0]))
    CAGCi = np.genfromtxt('Results/datToPlotCAGC/%s.dat'%(CAGC[i]))
    plt.figure(figsize=(30.0,8.0))
    CAGCiSX = (CAGCi[:,0])
    ax = plt.gca()
    ax1l, = ax.plot(CAGCi[:,0][::5],CAGCi[:,1][::5],color='midnightblue',linestyle='none',marker='X',markeredgecolor='midnightblue',markeredgewidth=None,markersize=5,alpha=1.0,label= '"%s"'%(CAGC[i]), zorder=3)
    ax2l, = ax.plot(CAGCi[:,0][::5],CAGCi[:,3][::5],color='darkorange',linestyle='none',marker='P',markeredgecolor='darkorange',markeredgewidth=None,markersize=5,alpha=1.0,label= '"%s"'%(CAGC[i]), zorder=4)
    ax3l, = ax.plot(CAGCi[:,0][::5],CAGCi[:,2][::5],color='c',linestyle='none',marker=6,markeredgecolor='c',markeredgewidth=None,markersize=5,alpha=1.0,label= '"%s"'%(CAGC[i]), zorder=5)
    ax4l, = ax.plot(CAGCi[:,0][::5],CAGCi[:,5][::5],color='m',linestyle='none',marker=7,markeredgecolor='m',markeredgewidth=None,markersize=5,alpha=1.0,label= '"%s"'%(CAGC[i]), zorder=6)
    ax.set_xlim([0,np.max(CAGCiSX)])
    ax.set_ylim([0.0,30.0])
    ax.yaxis.tick_left()
    ax.yaxis.set_label_position("left")
    ax.tick_params(axis='y',which='major',right=1,labelleft='on',labelright='on')
    plt.rcParams['xtick.labelsize'] = 30.0
    plt.rcParams['ytick.labelsize'] = 30.0
    plt.rc('font',weight='bold')
    ax.set_xlabel('Time [ns]',size=30.0)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    #for label in ax.xaxis.get_ticklabels()[1::2]:
     # label.set_visible(False)
    ax.xaxis.grid(True,'major',linewidth=2.0,color='k',zorder=0)
    ax.xaxis.grid(True,'minor',linewidth=1.0,color='gray',zorder=1)
    ax.yaxis.grid(True,'major',linewidth=2.0,color='k',zorder=0)
    ax.yaxis.grid(True,'minor',linewidth=1.0,color='gray',zorder=1)
    l = plt.legend([ax2l,ax1l,ax4l,ax3l],[r'$r_{(SC,SC)}[\mathrm{\AA}]\mid$'+r'$\langle x \rangle$='+'%.1f'%(CAGCstat[i][5])+r'$[\mathrm{\AA}]\mid$'+r' $SD_{x_{i}}$='+'%.1f'%(CAGCstat[i][6])+r'$[\mathrm{\AA}]\mid$', r'$r_{(CA,CA)}[\mathrm{\AA}]\mid$'+r'$\langle x \rangle$='+'%.1f'%(CAGCstat[i][1])+r'$[\mathrm{\AA}]\mid$'+r' $SD_{x_{i}}$='+'%.1f'%(CAGCstat[i][2])+r'$[\mathrm{\AA}]\mid$',r'$r_{(SC,MC)}[\mathrm{\AA}]\mid$'+r'$\langle x \rangle$='+'%.1f'%(CAGCstat[i][9])+r'$[\mathrm{\AA}]\mid$'+r' $SD_{x_{i}}$='+'%.1f'%(CAGCstat[i][10])+r'$[\mathrm{\AA}]$',r'$r_{(MC,MC)}[\mathrm{\AA}]\mid$'+r'$\langle x \rangle$='+'%.1f'%(CAGCstat[i][3])+r'$[\mathrm{\AA}]\mid$'+r' $SD_{x_{i}}$='+'%.1f'%(CAGCstat[i][4])+r'$[\mathrm{\AA}]$'],loc='best',numpoints=1,ncol=2,markerscale=3.0,prop={'size':30.0},borderpad=0.2,handletextpad=0.05)
    ax.text(-35,-4.5,'%s'%(CAGC[i]),size=(label_size),rotation=0,bbox=dict(facecolor='w',edgecolor='black'))
    #ax.text(700,-1.8,'Time 'r'$[\mathrm{ns}]$',size=label_size,fontsize=label_size)
    plt.savefig('Results/PlotCAGCvsTime/CAGCvsTime%s.svg'%(CAGC[i].replace(" ","")),format='svg',bbox_inches='tight',dpi=None,facecolor='#a9f971')
    plt.close()

print ('Finished CA/GC plots')

