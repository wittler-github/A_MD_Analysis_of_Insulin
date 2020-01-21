import matplotlib
import numpy as np
import re
from operator import itemgetter
import subprocess
print ("MPL version ")
print (matplotlib.__version__)


##
## SORT & COUNT HYDROGEN BONDS ##
##
print ('Commencing HB plots counting')

subprocess.call(["rm","-rf","Results/SortCountHB"]) 
subprocess.call(["mkdir","Results/SortCountHB"]) 

HB = np.genfromtxt('Results/datToPlotHB/HBperc_all_filenames.dat',dtype=('U50',None,None,None,None),delimiter="\t");

HBsortlist = [ [] for _ in range(len(HB))]

###CHAIN A B###

AchainHBlist = [ [] for _ in range(len(HB))]
BchainHBlist = [ [] for _ in range(len(HB))]
ABchainHBlist = [ [] for _ in range(len(HB))]

for i in range(len(HB)):
  HBresIDs = re.findall('[A]\(\d+|[D]\(\d+',HB[i][0])
  HBresIDs = ''.join(HBresIDs)
  HBresIDs = [int(s) for s in re.findall('\d+',HBresIDs)]
  HBsortlist[i].extend(HB[i])
  HBsortlist[i].extend(HBresIDs)
  print (HBsortlist[i])
  print ("\n")

HBsortlist = sorted(HBsortlist,key=itemgetter(5,6))

for i in range(len(HBsortlist)):
  print (HBsortlist[i])
  print ("HBsortlist\n")

for i in range(len(HBsortlist)):
  if ((HBsortlist[i][5]<22)&(HBsortlist[i][6]<22)):
   if not (re.search('[A-Z]\(\d+ [A-Z]{3} C',HBsortlist[i][0])):
     if (HBsortlist[i][2]>0.0):
       AchainHBlist[i].extend(HBsortlist[i])
AchainHBlist=[x for x in AchainHBlist if x != []]
for i in range(len(AchainHBlist)):
  print (AchainHBlist[i],'AchainHBlist\n')

for i in range(len(HBsortlist)):
  if ((HBsortlist[i][5]<52)&(HBsortlist[i][6]<52)&(HBsortlist[i][5]>21)&(HBsortlist[i][6]>21)):
   if not (re.search('[A-Z]\(\d+ [A-Z]{3} C',HBsortlist[i][0])):
     if (HBsortlist[i][2]>0.0):
       BchainHBlist[i].extend(HBsortlist[i])
BchainHBlist=[x for x in BchainHBlist if x != []]
for i in range(len(BchainHBlist)):
  print (BchainHBlist[i],'BchainHBlist\n')
    
for i in range(len(HBsortlist)):
  if (((HBsortlist[i][5]<52)&(HBsortlist[i][6]<52)&(HBsortlist[i][5]>21)&(HBsortlist[i][6]<22)))|((HBsortlist[i][5]<52)&(HBsortlist[i][6]<52)&((HBsortlist[i][5]<22)&(HBsortlist[i][6]>21))):
   if (not (re.search('[A-Z]\(\d+ [A-Z]{3} C',HBsortlist[i][0]))):
     if (HBsortlist[i][2]>0.0):
       ABchainHBlist[i].extend(HBsortlist[i])
ABchainHBlist=[x for x in ABchainHBlist if x != []]
for i in range(len(ABchainHBlist)):
  print (ABchainHBlist[i],'ABchainHBlist\n')

ABABchainHBlist=[]
ABABchainHBlist.extend(AchainHBlist)
ABABchainHBlist.extend(BchainHBlist)
ABABchainHBlist.extend(ABchainHBlist)

ABABchainHBlist=np.array(ABABchainHBlist,dtype=object)


hbcPerc = np.arange(0.0,5.0,5.0)
hbc = np.zeros(len(hbcPerc),dtype='U50,int,int,int,int,int') 
io=0
for i in range(len(hbcPerc)):
  hbc[io][0] = 'nrHBlt'+'%s'%(hbcPerc[i])+'p_A/B/AB(A+B+AB)' 
  hbc[io][1]=hbcPerc[i]
  hbc[io][2]=sum(x[2] >= hbcPerc[i] for x in AchainHBlist)
  hbc[io][3]=sum(x[2] >= hbcPerc[i] for x in BchainHBlist)
  hbc[io][4]=sum(x[2] >= hbcPerc[i] for x in ABchainHBlist)
  hbc[io][5]=sum(x[2] >= hbcPerc[i] for x in ABABchainHBlist)
  io+=1
  
  htpc_ABABchainHBlist = [x for x in ABABchainHBlist if x[2] >= hbcPerc[i]]
  htpc_lenHB=len(htpc_ABABchainHBlist)
  if htpc_lenHB !=0:
    np.savetxt('Results/SortCountHB/Sorted_chainaAB_HB_ht_'+'%d'%(hbcPerc[i])+'%'+'_nrHBs_%d'%(htpc_lenHB)+'.dat',htpc_ABABchainHBlist,fmt="%50s\t%-+.2f\t%-+.2f\t%-+.2f\t%-+.2f\t%d\t%d")

np.savetxt('Results/SortCountHB/Counted_chainAB_HB.dat',hbc,fmt="%40s\t%d\t%d\t%d\t%d\t%d\n")  

###CHAIN C D###

CchainHBlist = [ [] for _ in range(len(HB))]
DchainHBlist = [ [] for _ in range(len(HB))]
CDchainHBlist = [ [] for _ in range(len(HB))]


for i in range(len(HBsortlist)):
  if ((HBsortlist[i][5]>51)&(HBsortlist[i][6]>51)&(HBsortlist[i][5]<73)&(HBsortlist[i][6]<73)):
   if not (re.search('[A-Z]\(\d+ [A-Z]{3} C',HBsortlist[i][0])):
     if (HBsortlist[i][2]>0.0):
       CchainHBlist[i].extend(HBsortlist[i])
CchainHBlist=[x for x in CchainHBlist if x != []]
for i in range(len(CchainHBlist)):
  print (CchainHBlist[i],'CchainHBlist\n')

for i in range(len(HBsortlist)):
  if ((HBsortlist[i][5]>72)&(HBsortlist[i][6]>72)):
   if not (re.search('[A-Z]\(\d+ [A-Z]{3} C',HBsortlist[i][0])):
     if (HBsortlist[i][2]>0.0):
       DchainHBlist[i].extend(HBsortlist[i])
DchainHBlist=[x for x in DchainHBlist if x != []]
for i in range(len(DchainHBlist)):
  print (DchainHBlist[i],'DchainHBlist\n')
    
for i in range(len(HBsortlist)):
  if (((HBsortlist[i][5]>51)&(HBsortlist[i][6]>51)&(HBsortlist[i][5]>72)&(HBsortlist[i][6]<73)))|(((HBsortlist[i][5]>51)&(HBsortlist[i][6]>51)&(HBsortlist[i][5]<73)&(HBsortlist[i][6]>72))):
   if (not (re.search('[A-Z]\(\d+ [A-Z]{3} C',HBsortlist[i][0]))):
     if (HBsortlist[i][2]>0.0):
       CDchainHBlist[i].extend(HBsortlist[i])
CDchainHBlist=[x for x in CDchainHBlist if x != []]
for i in range(len(CDchainHBlist)):
  print (CDchainHBlist[i],'CDchainHBlist\n')

CDCDchainHBlist=[]
CDCDchainHBlist.extend(CchainHBlist)
CDCDchainHBlist.extend(DchainHBlist)
CDCDchainHBlist.extend(CDchainHBlist)

CDCDchainHBlist=np.array(CDCDchainHBlist,dtype=object)


hbcPerc = np.arange(0.0,5.0,5.0)
hbc = np.zeros(len(hbcPerc),dtype='U50,int,int,int,int,int') 
io=0
for i in range(len(hbcPerc)):
  hbc[io][0] = 'nrHBlt'+'%s'%(hbcPerc[i])+'p_C/D/CD(C+D+CD)' 
  hbc[io][1]=hbcPerc[i]
  hbc[io][2]=sum(x[2] >= hbcPerc[i] for x in CchainHBlist)
  hbc[io][3]=sum(x[2] >= hbcPerc[i] for x in DchainHBlist)
  hbc[io][4]=sum(x[2] >= hbcPerc[i] for x in CDchainHBlist)
  hbc[io][5]=sum(x[2] >= hbcPerc[i] for x in CDCDchainHBlist)
  io+=1
  
  htpc_CDCDchainHBlist = [x for x in CDCDchainHBlist if x[2] >= hbcPerc[i]]
  htpc_lenHB=len(htpc_CDCDchainHBlist)
  if htpc_lenHB !=0:
    np.savetxt('Results/SortCountHB/Sorted_chainCD_HB_ht_'+'%d'%(hbcPerc[i])+'%'+'_nrHBs_%d'%(htpc_lenHB)+'.dat',htpc_CDCDchainHBlist,fmt="%50s\t%-+.2f\t%-+.2f\t%-+.2f\t%-+.2f\t%d\t%d")

np.savetxt('Results/SortCountHB/Counted_chainCD_HB.dat',hbc,fmt="%40s\t%d\t%d\t%d\t%d\t%d\n")  


###CHAIN ABCD###

ABCDchainHBlist = [ [] for _ in range(len(HB))]

for i in range(len(HBsortlist)):
  if (((HBsortlist[i][5]>51)&(HBsortlist[i][6]<52))|((HBsortlist[i][5]<52)&(HBsortlist[i][6]>51))):
   if (not (re.search('[A-Z]\(\d+ [A-Z]{3} C',HBsortlist[i][0]))):
     if (HBsortlist[i][2]>0.0):
       ABCDchainHBlist[i].extend(HBsortlist[i])
ABCDchainHBlist=[x for x in ABCDchainHBlist if x != []]
for i in range(len(ABCDchainHBlist)):
  print (ABCDchainHBlist[i],'ABCDchainHBlist\n')

ABCDchainHBlist=np.array(ABCDchainHBlist,dtype=object)


hbcPerc = np.arange(0.0,5.0,5.0)
hbc = np.zeros(len(hbcPerc),dtype='U50,int,int') 
io=0
for i in range(len(hbcPerc)):
  hbc[io][0] = 'nrHBlt'+'%s'%(hbcPerc[i])+'p_ABCD' 
  hbc[io][1]=hbcPerc[i]
  hbc[io][2]=sum(x[2] >= hbcPerc[i] for x in ABCDchainHBlist)

  io+=1
  
  htpc_ABCDchainHBlist = [x for x in ABCDchainHBlist if x[2] >= hbcPerc[i]]
  htpc_lenHB=len(htpc_ABCDchainHBlist)
  if htpc_lenHB !=0:
    np.savetxt('Results/SortCountHB/Sorted_chainABCD_HB_ht_'+'%d'%(hbcPerc[i])+'%'+'_nrHBs_%d'%(htpc_lenHB)+'.dat',htpc_ABCDchainHBlist,fmt="%50s\t%-+.2f\t%-+.2f\t%-+.2f\t%-+.2f\t%d\t%d")

np.savetxt('Results/SortCountHB/Counted_chainABCD_HB.dat',hbc,fmt="%40s\t%d\t%d\n")  


###CHAIN ALL HB###

AllHBHBlist = [ [] for _ in range(len(HB))]

for i in range(len(HBsortlist)):
  if (1):
   if (not (re.search('[A-Z]\(\d+ [A-Z]{3} C',HBsortlist[i][0]))):
     if (HBsortlist[i][2]>0.0):
       AllHBHBlist[i].extend(HBsortlist[i])
AllHBHBlist=[x for x in AllHBHBlist if x != []]
for i in range(len(AllHBHBlist)):
  print (AllHBHBlist[i],'AllHBHBlist\n')

AllHBHBlist=np.array(AllHBHBlist,dtype=object)


hbcPerc = np.arange(0.0,5.0,5.0)
hbc = np.zeros(len(hbcPerc),dtype='U50,int,int') 
io=0
for i in range(len(hbcPerc)):
  hbc[io][0] = 'nrHBlt'+'%s'%(hbcPerc[i])+'p_ABCD' 
  hbc[io][1]=hbcPerc[i]
  hbc[io][2]=sum(x[2] >= hbcPerc[i] for x in AllHBHBlist)

  io+=1
  
  htpc_AllHBHBlist = [x for x in AllHBHBlist if x[2] >= hbcPerc[i]]
  htpc_lenHB=len(htpc_AllHBHBlist)
  if htpc_lenHB !=0:
    np.savetxt('Results/SortCountHB/Sorted_chain_AllHB_ht_'+'%d'%(hbcPerc[i])+'%'+'_nrHBs_%d'%(htpc_lenHB)+'.dat',htpc_AllHBHBlist,fmt="%50s\t%-+.2f\t%-+.2f\t%-+.2f\t%-+.2f\t%d\t%d")

np.savetxt('Results/SortCountHB/Counted_AllHB.dat',hbc,fmt="%40s\t%d\t%d\n")  


print ('END OF SCRIPT')
## END OF SCRIPT ##
