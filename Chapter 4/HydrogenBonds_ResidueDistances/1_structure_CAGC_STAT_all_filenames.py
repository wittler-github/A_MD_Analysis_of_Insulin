
import numpy as np

QFC = np.genfromtxt('Results/datToPlotCAGC/CAGC_all_filenames.dat',dtype='U20',delimiter="\n");

eqt = 0.0

meanQFCarray = np.zeros((7,len(QFC)))

ab = np.zeros(len(QFC),dtype=[('var1','U50'),('var2',float),('var3',float),('var4',float),('var5',float),('var6',float),('var7',float)]) 
ab['var1'] = QFC

for i in range(len(QFC)):
  QFCi = np.genfromtxt('Results/datToPlotCAGC/%s.dat'%(QFC[i]))
  nrXrow = 1
  nrYcol = len(QFCi)-1

  for Ycol in range(1,nrYcol+1):
    meanQFC = 0

    for j in range(np.int_(eqt),nrXrow):
      meanQFC = meanQFC+QFCi[Ycol]/(np.float(nrXrow)-eqt)  
    meanQFCarray[Ycol-1][i] = meanQFC
 
ab['var2'] = meanQFCarray[0]
ab['var3'] = meanQFCarray[1]
ab['var4'] = meanQFCarray[2]
ab['var5'] = meanQFCarray[3]
ab['var6'] = meanQFCarray[4]
ab['var7'] = meanQFCarray[5]

np.savetxt('Results/datToPlotCAGC/QFCstat_all_filenames.dat',ab,fmt="%20s\t%-+.6f\t%-+.6f\t%-+.6f\t%-+.6f\t%-+.6f\t%-+.6f")  

print ("END OF SCRIPT")
