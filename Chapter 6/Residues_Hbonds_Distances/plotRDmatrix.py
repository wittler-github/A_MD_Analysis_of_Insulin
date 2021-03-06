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

label_size1=11
label_size2=10

RDdata = np.genfromtxt('Results/datToPlotCAGC/QFCstat_all_filenames.dat',dtype=('U50',None,None,None,None,None,None,None,None,None,None),delimiter="\t");

ml = 51
RDM = (ml,ml)

RDM = np.zeros(RDM)

print (RDM)
print (RDdata)

for i in range(len(RDdata)):
  RDij = re.findall('(\d+)\s+[A-Z]{3}\s+to\s+(\d+)\s+[A-Z]{3}',RDdata[i][0])

  RDM[(int(RDij[0][0])-1)][(int(RDij[0][1])-1)] =  RDdata[i][3]
  RDM[(int(RDij[0][1])-1)][(int(RDij[0][0])-1)] =  RDdata[i][9]

  print (RDM[(int(RDij[0][0])-1)][(int(RDij[0][1])-1)])
  print (RDij[0][0]+' '+RDij[0][1])
  print ('I=%d'%(i)+'\n')
  
RDM_a = np.array(RDM)
Sum_RDM = np.sum(RDM_a)
print ('Sum_RDM = ',Sum_RDM)
print ('Maxij_RDM = ',np.max(RDM_a))
print ('Resij of Maxij =',np.argwhere(RDM_a.max() == RDM_a))

'''
### Plot the RDM less than 30 Å ###
fig = plt.figure()
ax = fig.gca()

cmap=plt.cm.rainbow_r
norm = matplotlib.colors.BoundaryNorm(np.arange(-1,31+1,1), cmap.N)
plt.imshow(RDM_a,cmap=cmap,norm=norm,extent=(0.5,51.5,51.5,0.5))
cbar = plt.colorbar(ticks=[0,5,10,15,20,25,30],drawedges=1,pad=0.01)
cbar.ax.tick_params(labelsize=label_size2) 

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
color_vrc = ['o','k','k']

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

plt.xlabel('A-chain/B-chain',size=label_size1)
plt.ylabel('A-chain/B-chain',size=label_size1)
plt.savefig('Results/PlotCAGC/RDM_SCSC_CACA_c1a5_lt30Å.svg',dpi=300,bbox_inches='tight',pad_inches=0.01)
plt.clf()

plt.close('all')
'''
### Plot the RDM less than 10 ###
fig = plt.figure()
ax = fig.gca()

cmap=plt.cm.cubehelix
norm = matplotlib.colors.BoundaryNorm(np.arange(-0.5,10+1,0.5), cmap.N)
plt.imshow(RDM_a,cmap=cmap,norm=norm,extent=(0.5,51.5,51.5,0.5))
cbar = plt.colorbar(ticks=np.linspace(0,10,20+1),drawedges=1,pad=0.01)
cbar.ax.tick_params(labelsize=label_size2) 

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

plt.xlabel('A-chain/B-chain',size=label_size1)
plt.ylabel('A-chain/B-chain',size=label_size1)
plt.savefig('Results/PlotCAGC/RDM_SCMC_MCMC_c3a9_lt10Å.svg',dpi=300,bbox_inches='tight',pad_inches=0.01)
plt.clf()

plt.close('all')

print("END OF SCRIPT")

