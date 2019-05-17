# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 14:10:28 2016

@author: se
"""

#import of modules

import csv
import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import griddata

#import of Gamry .DTA file with different paragraphs of CV data and a long QCM data paragraph in the end

u = '/home/se/Dropbox/Masterarbeit/Projekt/data/raw/eQCM/C_beschichtet/#2/2016-05-19/QCM_RCVv8_10mV_all/QCM_RCVv8_10mV_all.DTA' #location of .DTA file
f = open(u, 'r') #open .DTA file
origin = os.path.dirname(u)+'/' #folder of .DTA file
data=[] #initiate the list/tuple

for line in f:
    line = line.strip()
    columns = line.split()
    data.append(columns)

data = filter(None, data) #removes empty lines in .dta file
    
cvlist = [] #starting points of CV cycles
qcmlist = [] #starting point of QCM data 
qcmstart = [] #time the QCM was measuring before the electrochemical measurements were started and postition

for i in data:
    if i[0].startswith('CURVE'):
            #print data.index(i)
            cvlist.append(data.index(i))
           
    if  'QCMOFFSETTIME'in i:
            qcmstart.append(data[data.index(i)][2])
            qcmstart.append(data.index(i)) #time the QCM was measuring before the electrochemical measurements were started in sec.
    if  'QCMCurve' in i:
            qcmlist = data.index(i)
            #print 'qcm position found:'           
            #print qcmlist
            #print data.index(i)    
            
qcm_time_zero = float((qcmstart[0]).replace(',','.'))                        
#print 'finished'
print ('Es wurden '+repr(len(cvlist)) +' CV Zyklen gefunden.')
print ('Diese starten in den Zeilen: '+ repr(cvlist))
print ('Die QCM Daten beginnen ab Zeile: '+ repr(qcmlist))

header = data[0:(cvlist[0])]

#CV Zyklen splitten:

cvname = [] #initiate cv data names

a=0 
for i in cvlist:
    a = a+1
    cvname.append('CV_{0:02d}'.format(a))
    
cvdata = []    #initiate echem data collection, enthält alle CV Zyklen(getrennt)

for i in cvlist:
    if cvlist.index(i) < len(cvlist)-1:
        cvdata.append(data[cvlist[cvlist.index(i)]+1:cvlist[cvlist.index(i)+1]])
        cvdata[cvlist.index(i)][1][3] = 'V_vs._Ref.'
        cvdata[cvlist.index(i)][1][4:6] = ''
    else:
        cvdata.append(data[cvlist[cvlist.index(i)]+1:(qcmstart[1]-2)])
        cvdata[cvlist.index(i)][1][3] = 'V_vs._Ref.'
        cvdata[cvlist.index(i)][1][4:6] = ''
        
for i in range(len(cvdata)):
    for j in range(len(cvdata[i])):
        for k in range(len(cvdata[i][j])-1):
           cvdata[i][j][k] = cvdata[i][j][k].replace(',','.')
        
qcmdata = data[qcmlist+1:len(data)]

for i in range(len(qcmdata)):
    for j in range(len(qcmdata[i])):
          qcmdata[i][j] = qcmdata[i][j].replace(',','.')

#saving QCM and CV data to files
nr = 0

#save header

nameheader = origin + 'header.txt'
file = open(nameheader, 'w')
wr = csv.writer(file, quoting=csv.QUOTE_ALL)
wr.writerows((header))
file.close()


for i in range(len(cvname)):
    namecv = origin + cvname[i] + '.txt'
    cv_loop =np.asarray(cvdata[i][2:])
    cv_loop = np.asarray(cv_loop[:][:,1:5], dtype = np.float32)
    np.savetxt(namecv, cv_loop, delimiter=',', header='T/s, Q/C, U_ref/V, I/A')
    #file = open(namecv, 'wb')
    #wr = csv.writer(file, quoting=csv.QUOTE_ALL)
    #wr.writerows((cvdata[nr]))
    #file.close()
    #nr = nr +1

cvtime = [] #strating time of each cv cycle
for i in range(len(cvname)):
    cvtime.append(float(cvdata[i][2][1]))



#modify QCM data to start at 0 seconds    
start = 0 

while start < len(qcmdata):
      if start <2:
          start= start + 1
      else:
          qcmdata[start][1] = float(qcmdata[start][1]) - qcm_time_zero
          start = start + 1
     
   
nameqcm = origin + 'qcm_data.txt'
file = open(nameqcm, 'w')
for i in range(len(qcmdata)):
  file.write("%s \n" % qcmdata[i])
file.close()

# End of wrtining original data to seperate files

q = 0
qcmname = []
for i in cvname:
   q = q+1
   qcmname.append('QCM_{0:02d}'.format(q))


eqcmdata = np.asarray(qcmdata)[2:][:,1:4]
eqcmdata = np.asarray(eqcmdata, dtype=np.float32)
eqcmname = origin + 'eqcm_data.txt'
np.savetxt(eqcmname, eqcmdata, delimiter=',', header='T/s , fs/MHz , fp/MHz')

#split qcmdata to CV cycles
t = 0
qcmtime = []
for i in cvtime:
  if eqcmdata[:][t][0] == 0.0:
    qcmtime.append(t)
    t = t+1
  else:
    while eqcmdata[:][t][0] < i:      
        t = t+1
    qcmtime.append(t)
    t = t+1

l = 0
for i in range(len(qcmtime)):
    if qcmtime[i]<qcmtime[-1]:
        np.savetxt(origin + qcmname[i] + '.txt', eqcmdata[l:qcmtime[i+1]], delimiter=',', header='T/s , fs/MHz , fp/MHz')
        l = qcmtime[i+1]
    else:
        np.savetxt(origin +qcmname[i] + '.txt', eqcmdata[qcmtime[-1]:], delimiter=',', header='T/s , fs/MHz , fp/MHz')
        
  
for i in range(len(cvname)):
    CV = np.loadtxt(origin + cvname[i] + '.txt', delimiter=',', skiprows=1)
    QCM = np.loadtxt(origin + qcmname[i] + '.txt', delimiter=',', skiprows=2, usecols =(1,))
    #QCMtime = QCMtime.flatten
    time = CV[:,[0]]
    time = time.flatten()
    Q = CV[:,[1]]
    Q = Q.flatten()
    U = CV[:,[2]]
    U = U.flatten()
    I = CV[:,[3]]
    I = I.flatten()
    plt.plot(time, U)
    plt.ylabel('U / V')
    plt.xlabel('time / s')
    savefig(origin + 'U_vs_t_{0:02d}.pdf'.format(i+1), bbox_inches='tight')
    plt.clf()
    plt.close()   
    fitQ = interp1d(time, Q)
    fitU = interp1d(time, U)
    fitI = interp1d(time, I)
    
    
    
    
    
    
    
    
    