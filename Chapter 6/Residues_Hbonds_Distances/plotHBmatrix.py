import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc
import matplotlib as mpl
import re
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
import math
#print ("MPL version ")
#print (matplotlib.__version__)

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode']=True

fontsize_1=10
label_size1=10
label_size2=12

HBdata = np.genfromtxt('Results/SortCountHB/Sorted_HB_ht_25%_nrHBs_28.dat',dtype=('U50',None,None,None,None,None,None),delimiter="\t");

ml = 51
HBM = (ml,ml)

HBM = np.zeros(HBM)

print (HBM)
print (HBdata)

for i in range(len(HBdata)):
  HBij = re.findall('[D|A]\((\d+)',HBdata[i][0])
  N_HBij = re.findall('[D|A]\(\d+\s+[A-Z]+\s+([A-Z]+)',HBdata[i][0])

  if ((N_HBij[0]!='N')&(N_HBij[0]!='C')&(N_HBij[0]!='O')&(N_HBij[0]!='OT1')&(N_HBij[0]!='OT2'))|((N_HBij[1]!='N')&(N_HBij[1]!='C')&(N_HBij[1]!='O')&(N_HBij[1]!='OT1')&(N_HBij[1]!='OT2'))&((int(HBij[0])-1)!=(int(HBij[1])-1)):
    if (int(HBij[1])>int(HBij[0])):
      HBM[(int(HBij[1])-1)][(int(HBij[0])-1)] +=  1  
    if (int(HBij[1])<int(HBij[0])):
      HBM[(int(HBij[0])-1)][(int(HBij[1])-1)] +=  1   

  if ((N_HBij[0]=='N')|(N_HBij[0]=='C')|(N_HBij[0]=='O')|(N_HBij[0]=='OT1')|(N_HBij[0]=='OT2'))&((N_HBij[1]=='N')|(N_HBij[1]=='C')|(N_HBij[1]=='O')|(N_HBij[1]=='OT1')|(N_HBij[1]=='OT2'))&((int(HBij[0])-1)!=(int(HBij[1])-1)):
    if (int(HBij[0])>int(HBij[1])):
      HBM[(int(HBij[1])-1)][(int(HBij[0])-1)] +=  1  
    if (int(HBij[0])<int(HBij[1])):
      HBM[(int(HBij[0])-1)][(int(HBij[1])-1)] +=  1   
   
  if ((int(HBij[0])-1)==(int(HBij[1])-1)):
    HBM[(int(HBij[0])-1)][(int(HBij[1])-1)] +=  1  

  print (HBM[(int(HBij[0])-1)][(int(HBij[1])-1)])
  print (HBij[0]+' '+HBij[1])
  print ('I=%d'%(i)+'\n')
  
HBM_a = np.array(HBM)
Sum_HBM = np.sum(HBM_a)
print ('Sum_HBM = ',Sum_HBM)
print ('Maxij_HBM = ',np.max(HBM_a))
print ('Resij of Maxij =',np.argwhere(HBM_a.max() == HBM_a))


### Plot the HBM less than 10 ###
fig = plt.figure()
ax = fig.gca()

#cmap=cmap=mcolors.ListedColormap(['w','blue'])
#cmap=cmap=mcolors.ListedColormap(['w','blue','red'])
#cmap=cmap=mcolors.ListedColormap(['w','blue','red','limegreen'])
cmap=cmap=mcolors.ListedColormap(['w','blue','red','limegreen','slategrey'])
#cmap=cmap=mcolors.ListedColormap(['w','blue','red','limegreen','slategrey','goldenrod'])
#cmap=cmap=mcolors.ListedColormap(['w','blue','red','limegreen','slategrey','goldenrod','hotpink'])
#cmap=cmap=mcolors.ListedColormap(['w','blue','red','limegreen','slategrey','goldenrod','hotpink','olive','darkorchid','k','c'])
#cmap=cmap=mcolors.ListedColormap(['w','blue','red','limegreen','slategrey','goldenrod','hotpink','olive','darkorchid','k','c','yellow'])
norm = matplotlib.colors.BoundaryNorm(np.arange(-0.5,math.ceil(np.max(HBM_a))+1,1),cmap.N)
plt.imshow(HBM_a,cmap=cmap,norm=norm,extent=(0.5,51+0.5,51+0.5,0.5))
cbar = plt.colorbar(ticks=np.linspace(0,math.ceil(np.max(HBM_a)),math.ceil(np.max(HBM_a))+1),drawedges=1,pad=0.01)
cbar.ax.tick_params(labelsize=label_size1) 

ax.invert_yaxis()
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
EI_HBx = ['(1 GLY 1)','(2 ILE 2)','(3 VAL 3)','(4 GLU 4)','(5 GLN 5)','(6 CYS 6)','(7 CYS 7)','(8 THR 8)','(9 SER 9)','(10 ILE 10)','(11 CYS 11)','(12 SER 12)','(13 LEU 13)','(14 TYR 14)','(15 GLN 15)','(16 LEU 16)','(17 GLU 17)','(18 ASN 18)','(19 TYR 19)','(20 CYS 20)','(21 ASN 21)','(22 PHE 1)','(23 VAL 2)','(24 ASN 3)','(25 GLN 4)','(26 HIS 5)','(27 LEU 6)','(28 CYS 7)','(29 GLY 8)','(30 SER 9)','(31 HIS 10)','(32 LEU 11)','(33 VAL 12)','(34 GLU 13)','(35 ALA 14)','(36 LEU 15)','(37 TYR 16)','(38 LEU 17)','(39 VAL 18)','(40 CYS 19)','(41 GLY 20)','(42 GLU 21)','(43 ARG 22)','(44 GLY 23)','(45 PHE 24)','(46 PHE 25)','(47 TYR 26)','(48 THR 27)','(49 LYS 28)','(50 PRO 29)','(51 THR 30)']
EI_HBx=np.append(EI_HBx[0],EI_HBx)
ax.set_xticklabels(EI_HBx)
ax.set_yticklabels(EI_HBx)
ax.tick_params(axis='x',which='major',labelsize=5,pad=0.05)
ax.tick_params(axis='y',which='major',labelsize=5,pad=-2.7)
for tick in ax.get_xticklabels():
  tick.set_rotation(-90)

vrC = [0,21,51]
color_vrc = ['o','gold','lawngreen']

for i in range(1,len(vrC)):
  k=1
  for j in range(vrC[i-1],vrC[i]):
    print (j)
    if (k==1.0)|(j==(vrC[i]-1)):
      ax.axhline((j+1), linestyle='-', color=color_vrc[i],linewidth=1.0)
      ax.axvline((j+1), linestyle='-', color=color_vrc[i],linewidth=1.0)
    if (0.0==((j+1)%5)):
      ax.axhline((j+1), linestyle='-.', color=color_vrc[i],linewidth=0.5)
      ax.axvline((j+1), linestyle='-.', color=color_vrc[i],linewidth=0.5)
    if (0.0==((j+1)%10)):
      ax.axhline((j+1), linestyle='-', color=color_vrc[i],linewidth=0.5)
      ax.axvline((j+1), linestyle='-', color=color_vrc[i],linewidth=0.5)
    if (k!=1.0)|(j!=(vrC[i]-1)):
      ax.axhline((j+1), linestyle=':', color=color_vrc[i],linewidth=0.5)
      ax.axvline((j+1), linestyle=':', color=color_vrc[i],linewidth=0.5)
    k += 1


plt.savefig('Results/PlotHB/HBM_SCSCaSCMC_MCMC_vplt30_(Sorted_HB_ht_25%_nrHBs_28)_dpi300.svg',dpi=300,bbox_inches='tight',pad_inches=0.01)
plt.clf()

plt.close('all')

print("END OF SCRIPT")

