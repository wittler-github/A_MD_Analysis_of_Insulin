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
import matplotlib.gridspec as gridspec
#print ("MPL version ")
#print (matplotlib.__version__)

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode']=True
mpl.rcParams["figure.figsize"] = [6.4,6.4]


HBdata = np.genfromtxt('Results/SortCountHB/Sorted_chain_AllHB_ht_0%_nrHBs_190.dat',dtype=('U50',None,None,None,None,None,None),delimiter="\t");

ml = 102
HBM = (ml,ml)

HBM = np.zeros(HBM)

print (HBM)
print (HBdata)

miv0 = 1
miv1 = 1

for i in range(len(HBdata)):
  HBij = re.findall('[D|A]\((\d+)',HBdata[i][0])
  N_HBij = re.findall('[D|A]\(\d+\s+[A-Z]+\s+([A-Z]+)',HBdata[i][0])

  if ((N_HBij[0]!='N')&(N_HBij[0]!='C')&(N_HBij[0]!='O')&(N_HBij[0]!='OT1')&(N_HBij[0]!='OT2'))|((N_HBij[1]!='N')&(N_HBij[1]!='C')&(N_HBij[1]!='O')&(N_HBij[1]!='OT1')&(N_HBij[1]!='OT2'))&((int(HBij[0])-1)!=(int(HBij[1])-1)):
    if (int(HBij[1])>int(HBij[0])):
      HBM[(int(HBij[1])-miv1)][(int(HBij[0])-miv0)] +=  1  
    if (int(HBij[1])<int(HBij[0])):
      HBM[(int(HBij[0])-miv0)][(int(HBij[1])-miv1)] +=  1   

  if ((N_HBij[0]=='N')|(N_HBij[0]=='C')|(N_HBij[0]=='O')|(N_HBij[0]=='OT1')|(N_HBij[0]=='OT2'))&((N_HBij[1]=='N')|(N_HBij[1]=='C')|(N_HBij[1]=='O')|(N_HBij[1]=='OT1')|(N_HBij[1]=='OT2'))&((int(HBij[0])-1)!=(int(HBij[1])-1)):
    if (int(HBij[0])>int(HBij[1])):
      HBM[(int(HBij[1])-miv1)][(int(HBij[0])-miv0)] +=  1  
    if (int(HBij[0])<int(HBij[1])):
      HBM[(int(HBij[0])-miv0)][(int(HBij[1])-miv1)] +=  1   
   
  if ((int(HBij[0])-1)==(int(HBij[1])-1)):
    HBM[(int(HBij[0])-miv0)][(int(HBij[1])-miv1)] +=  1  

  print (HBM[(int(HBij[0])-miv0)][(int(HBij[1])-miv1)])
  print (HBij[0]+' '+HBij[1])
  print ('I=%d'%(i)+'\n')

HBM_a = np.array(HBM)
Sum_HBM = np.sum(HBM_a)
print ('Sum_HBM = ',Sum_HBM)
print ('Maxij_HBM = ',np.max(HBM_a))
print ('Resij of Maxij =',np.argwhere(HBM_a.max() == HBM_a))


### Plot the count HBM ###

f,(ax1to4) = plt.subplots(2,2,sharex=True,sharey=True)
ax = f.gca()
f.subplots_adjust(wspace=0.027,hspace=0.027)

cmap=cmap=mcolors.ListedColormap(['w','blue','red'])
norm = matplotlib.colors.BoundaryNorm(np.arange(-0.5,math.ceil(np.max(HBM_a))+1,1),cmap.N)

im1 = ax1to4[0,0].imshow(HBM_a[51:102,0:51],extent=(0.5,51.5,51.5,0.5),cmap=cmap,norm=norm,aspect='equal')
#ax1to4[0,0].text(25, 25, 'UL')

im2 = ax1to4[1,0].imshow(HBM_a[0:51,0:51],extent=(0.5,51.5,51.5,0.5),cmap=cmap,norm=norm,aspect='equal')
#ax1to4[1,0].text(25, 25, 'LL')

im3 = ax1to4[0,1].imshow(HBM_a[51:102,51:102],extent=(0.5,51.5,51.5,0.5),cmap=cmap,norm=norm,aspect='equal')
#ax1to4[0,1].text(25, 25, 'UR')

im4 = ax1to4[1,1].imshow(HBM_a[0:51,51:102],extent=(0.5,51.5,51.5,0.5),cmap=cmap,norm=norm,aspect='equal')
#ax1to4[1,1].text(25, 25, 'LR')

cbar_ax = f.add_axes([0.91, 0.11, 0.033, 0.77])
cb = f.colorbar(im1,ticks=np.linspace(0,math.ceil(np.max(HBM_a)),math.ceil(np.max(HBM_a))+1),cax=cbar_ax,drawedges=1)

ax.invert_yaxis()
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
EI_HBx = ['(1 GLY 1)','(2 ILE 2)','(3 VAL 3)','(4 GLU 4)','(5 GLN 5)','(6 CYS 6)','(7 CYS 7)','(8 THR 8)','(9 SER 9)','(10 ILE 10)','(11 CYS 11)','(12 SER 12)','(13 LEU 13)','(14 TYR 14)','(15 GLN 15)','(16 LEU 16)','(17 GLU 17)','(18 ASN 18)','(19 TYR 19)','(20 CYS 20)','(21 ASN 21)','(22 PHE 1)','(23 VAL 2)','(24 ASN 3)','(25 GLN 4)','(26 HIS 5)','(27 LEU 6)','(28 CYS 7)','(29 GLY 8)','(30 SER 9)','(31 HIS 10)','(32 LEU 11)','(33 VAL 12)','(34 GLU 13)','(35 ALA 14)','(36 LEU 15)','(37 TYR 16)','(38 LEU 17)','(39 VAL 18)','(40 CYS 19)','(41 GLY 20)','(42 GLU 21)','(43 ARG 22)','(44 GLY 23)','(45 PHE 24)','(46 PHE 25)','(47 TYR 26)','(48 THR 27)','(49 PRO 28)','(50 LYS 29)','(51 ALA 30)']
EI_HBy = ['(1 GLY 1)','(2 ILE 2)','(3 VAL 3)','(4 GLU 4)','(5 GLN 5)','(6 CYS 6)','(7 CYS 7)','(8 THR 8)','(9 SER 9)','(10 ILE 10)','(11 CYS 11)','(12 SER 12)','(13 LEU 13)','(14 TYR 14)','(15 GLN 15)','(16 LEU 16)','(17 GLU 17)','(18 ASN 18)','(19 TYR 19)','(20 CYS 20)','(21 ASN 21)','(1 PHE 22)','(2 VAL 23)','(3 ASN 24)','(4 GLN 25)','(5 HIS 26)','(6 LEU 27)','(7 CYS 28)','(8 GLY 29)','(9 SER 30)','(10 HIS 31)','(11 LEU 32)','(12 VAL 33)','(13 GLU 34)','(14 ALA 35)','(15 LEU 36)','(16 TYR 37)','(17 LEU 38)','(18 VAL 39)','(19 CYS 40)','(20 GLY 41)','(21 GLU 42)','(22 ARG 43)','(23 GLY 44)','(24 PHE 45)','(25 PHE 46)','(26 TYR 47)','(27 THR 48)','(28 PRO 49)','(29 LYS 50)','(30 ALA 51)']
EI_HBx=np.append(EI_HBx[0],EI_HBx)
ax1to4[1,0].set_xticklabels(EI_HBx)
ax1to4[1,0].set_yticklabels(EI_HBx)
ax1to4[1,0].tick_params(axis='x',which='major',labelsize=3.25,pad=0.05)
for tick in ax1to4[1,0].get_xticklabels():
  tick.set_rotation(-90)
ax1to4[1,0].tick_params(axis='y',which='major',labelsize=3.25,pad=-2)
ax1to4[1,1].set_xticklabels(EI_HBx)
ax1to4[1,1].tick_params(axis='x',which='major',labelsize=3.25,pad=0.05)
for tick in ax1to4[1,1].get_xticklabels():
  tick.set_rotation(-90)
EI_HBy2=np.append(EI_HBy[0],EI_HBy)
ax1to4[0,0].set_yticklabels(EI_HBy2)
ax1to4[0,0].tick_params(axis='y',which='major',labelsize=3.25,pad=-2)
'''
for i in [0,1]:
  ax1to4[0,i].tick_params(axis='both',bottom=0,top=0,left=0,right=0)
  ax1to4[1,i].tick_params(axis='both',bottom=0,top=0,left=0,right=0)

ax1to4[1,0].tick_params(axis='x',bottom=1,length=2)
ax1to4[1,0].tick_params(axis='y',left=1,length=2)
ax1to4[1,1].tick_params(axis='x',bottom=1,length=2)
ax1to4[0,0].tick_params(axis='y',left=1,length=2)
'''
vrC = [0,21,51]
color_vrc1 = ['o','cornflowerblue','darkorange']
color_vrc2 = ['o','darkorchid','mediumseagreen']
for i in range(1,len(vrC)):
  k=1
  for j in range(vrC[i-1],vrC[i]):
    print (j)
    if (k==1.0)|(j==(vrC[i]-1)):
      ax1to4[0,0].axhline((j+1), linestyle='-', color=color_vrc1[i],linewidth=0.75)
      ax1to4[0,0].axvline((j+1), linestyle='-', color=color_vrc2[i],linewidth=0.75)
      ax1to4[0,1].axhline((j+1), linestyle='-', color=color_vrc1[i],linewidth=0.75)
      ax1to4[0,1].axvline((j+1), linestyle='-', color=color_vrc1[i],linewidth=0.75)
      ax1to4[1,0].axhline((j+1), linestyle='-', color=color_vrc2[i],linewidth=0.75)
      ax1to4[1,0].axvline((j+1), linestyle='-', color=color_vrc2[i],linewidth=0.75)
      ax1to4[1,1].axhline((j+1), linestyle='-', color=color_vrc2[i],linewidth=0.75)
      ax1to4[1,1].axvline((j+1), linestyle='-', color=color_vrc1[i],linewidth=0.75)
    if (0.0==((j+1)%5)):
      ax1to4[0,0].axhline((j+1), linestyle='-.', color=color_vrc1[i],linewidth=0.375)
      ax1to4[0,0].axvline((j+1), linestyle='-.', color=color_vrc2[i],linewidth=0.375)
      ax1to4[0,1].axhline((j+1), linestyle='-.', color=color_vrc1[i],linewidth=0.375)
      ax1to4[0,1].axvline((j+1), linestyle='-.', color=color_vrc1[i],linewidth=0.375)
      ax1to4[1,0].axhline((j+1), linestyle='-.', color=color_vrc2[i],linewidth=0.375)
      ax1to4[1,0].axvline((j+1), linestyle='-.', color=color_vrc2[i],linewidth=0.375)
      ax1to4[1,1].axhline((j+1), linestyle='-.', color=color_vrc2[i],linewidth=0.375)
      ax1to4[1,1].axvline((j+1), linestyle='-.', color=color_vrc1[i],linewidth=0.375)
    if (0.0==((j+1)%10)):
     ax1to4[0,0].axhline((j+1), linestyle='-', color=color_vrc1[i],linewidth=0.375)
     ax1to4[0,0].axvline((j+1), linestyle='-', color=color_vrc2[i],linewidth=0.375)
     ax1to4[0,1].axhline((j+1), linestyle='-', color=color_vrc1[i],linewidth=0.375)
     ax1to4[0,1].axvline((j+1), linestyle='-', color=color_vrc1[i],linewidth=0.375)
     ax1to4[1,0].axhline((j+1), linestyle='-', color=color_vrc2[i],linewidth=0.375)
     ax1to4[1,0].axvline((j+1), linestyle='-', color=color_vrc2[i],linewidth=0.375)
     ax1to4[1,1].axhline((j+1), linestyle='-', color=color_vrc2[i],linewidth=0.375)
     ax1to4[1,1].axvline((j+1), linestyle='-', color=color_vrc1[i],linewidth=0.375)
    if (k!=1.0)|(j!=(vrC[i]-1)):
      ax1to4[0,0].axhline((j+1), linestyle=':', color=color_vrc1[i],linewidth=0.375)
      ax1to4[0,0].axvline((j+1), linestyle=':', color=color_vrc2[i],linewidth=0.375)
      ax1to4[0,1].axhline((j+1), linestyle=':', color=color_vrc1[i],linewidth=0.375)
      ax1to4[0,1].axvline((j+1), linestyle=':', color=color_vrc1[i],linewidth=0.375)
      ax1to4[1,0].axhline((j+1), linestyle=':', color=color_vrc2[i],linewidth=0.375)
      ax1to4[1,0].axvline((j+1), linestyle=':', color=color_vrc2[i],linewidth=0.375)
      ax1to4[1,1].axhline((j+1), linestyle=':', color=color_vrc2[i],linewidth=0.375)
      ax1to4[1,1].axvline((j+1), linestyle=':', color=color_vrc1[i],linewidth=0.375)
    k += 1

ax1to4[1,0].set_xlabel('\ M2 AC/BC')
ax1to4[1,0].set_ylabel('\ M2 AC/BC')
ax1to4[1,1].set_xlabel('\ M1 AC/BC')
ax1to4[0,0].set_ylabel('\ M1 AC/BC')

plt.savefig('Results/PlotHB/HBM_Sorted_chain_AllHB_ht_0%_nrHBs_190.svg',dpi=300,bbox_inches='tight',pad_inches=0.025,facecolor='w')
plt.clf()

print("END OF SCRIPT")
