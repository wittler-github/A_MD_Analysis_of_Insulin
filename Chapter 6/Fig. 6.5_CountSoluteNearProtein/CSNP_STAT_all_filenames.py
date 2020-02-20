
import numpy as np

QFC = np.genfromtxt('Results/datToPlotCSNP/CSNP_all_filenames.dat',dtype=str);

eqt = 9.0

meanQFCarray = np.zeros((9,len(QFC)))
sdQFCarray = np.zeros((9,len(QFC)))

ab = np.zeros(len(QFC),dtype=[('var1','U50'),('var2',float),('var3',float),('var4',float),('var5',float),('var6',float),('var7',float),('var8',float),('var9',float)]) 
ab['var1'] = QFC

for i in range(len(QFC)):
  QFCi = np.genfromtxt('Results/datToPlotCSNP/%s.dat'%(QFC[i]))
  nrXrow = len(QFCi[:,0])
  nrYcol = len(QFCi[0,:])-1

  for Ycol in range(1,nrYcol+1):
    meanQFC = 0
    sdQFC = 0

    for j in range(np.int_(eqt),nrXrow):
      meanQFC = meanQFC+QFCi[j,Ycol]/(np.float(nrXrow)-eqt)  
    meanQFCarray[Ycol-1][i] = meanQFC
 
    for k in range(np.int_(eqt),nrXrow):
      sdQFC = sdQFC+np.power((QFCi[k,Ycol]-meanQFC),2)/((np.float(nrXrow)-eqt)-1.0) 
    sdQFC = np.sqrt(sdQFC)
    sdQFCarray[Ycol-1][i] = sdQFC

ab['var2'] = meanQFCarray[0]
ab['var3'] = sdQFCarray[0]
ab['var4'] = meanQFCarray[1]
ab['var5'] = sdQFCarray[1]
ab['var6'] = meanQFCarray[2]
ab['var7'] = sdQFCarray[2]
ab['var8'] = meanQFCarray[3]
ab['var9'] = sdQFCarray[3]

np.savetxt('Results/datToPlotCSNP/QFCstat_all_filenames.dat',ab,fmt="%12s\t%-+.6f\t%-+.6f\t%-+.6f\t%-+.6f\t%-+.6f\t%-+.6f\t%-+.6f\t%-+.6f")