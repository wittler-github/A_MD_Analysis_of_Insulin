
import numpy as np

QFC = np.genfromtxt('Results/datToPlotDA/E_DA_all_filenames.dat',dtype=('U10'),delimiter="\n");

eqt = 0.0

fracQFCarray = np.zeros((7,len(QFC)))

ab = np.zeros(len(QFC),dtype=[('var1','U10'),('var2',float),('var3',float),('var4',float),('var5',float),('var6',float),('var7',float),('var8',float)]) 
T=np.array(range(1,52))
ab['var1'] = T

for i in range(len(QFC)):
  QFCi = np.genfromtxt('Results/datToPlotDA/%s.dat'%(QFC[i]))
  nrXrow = 1
  nrYcol = len(QFCi)-1

  for Ycol in range(1,nrYcol+1):
    fracQFC = 0
    sdQFC = 0

    for j in range(np.int_(eqt),nrXrow):
      fracQFC = fracQFC+QFCi[Ycol]/(np.float(nrXrow)-eqt)  
    fracQFCarray[Ycol-1][i] = fracQFC

ab['var2'] = fracQFCarray[0]
ab['var3'] = fracQFCarray[1]
ab['var4'] = fracQFCarray[2]
ab['var5'] = fracQFCarray[3]
ab['var6'] = fracQFCarray[4]
ab['var7'] = fracQFCarray[5]
ab['var8'] = fracQFCarray[6]

np.savetxt('Results/datToPlotDA/E_QFC_DA_stat_all_filenames.dat',ab,fmt="%10s\t%-+.6f\t%-+.6f\t%-+.6f\t%-+.6f\t%-+.6f\t%-+.6f\t%-+.6f")  
