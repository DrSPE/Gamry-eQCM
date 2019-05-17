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
    cvname.append('CV'+repr(a))
    
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
    for j in range(len(cvdata[i])-1):
        for k in range(len(cvdata[i][j])-1):
           cvdata[i][j][k] = cvdata[i][j][k].replace(',','.')
        
qcmdata = data[qcmlist+1:len(data)]

for i in range(len(qcmdata)):
    for j in range(len(qcmdata[i])):
          qcmdata[i][j] = qcmdata[i][j].replace(',','.')

#saving QCM and CV data to files
nr = 0

#save header

nameheader = origin + 'header.csv'
file = open(nameheader, 'w')
wr = csv.writer(file, quoting=csv.QUOTE_ALL)
wr.writerows((header))
file.close()


for i in cvname:
    namecv = origin + cvname[cvname.index(i)] + '.csv'
    file = open(namecv, 'wb')
    wr = csv.writer(file, quoting=csv.QUOTE_ALL)
    wr.writerows((cvdata[nr]))
    file.close()
    nr = nr +1


#modify QCM data to start at 0 seconds    
start = 0 

while start < len(qcmdata):
      if start <2:
          start= start + 1
      else:
          qcmdata[start][1] = float(qcmdata[start][1]) - qcm_time_zero
          start = start + 1
     
     
nameqcm = origin + 'qcm_data.csv'
file = open(nameqcm, 'w')
wr = csv.writer(file, quoting=csv.QUOTE_ALL)
wr.writerows((qcmdata))
file.close()





#define fit function to align data to the same time axis 

#for i in eqcmdata['T'][data.index(i)]:
       #data['T_new] = (data['T'][data.index(i)]) - 17.09
       #return data['T_new]

    
df_qcm = pd.read_csv(nameqcm)
eqcmdata = df_qcm.iloc[:,1:4]

test = []
a=1
b=1
for i in cvname:
    namecvfile = origin + i + '.csv'
    df_cv0 = pd.read_csv(namecvfile)
    df_cv = df_cv0.iloc[:,1:5]
    cvheader = list(df_cv.columns.values)
    for j in range(len(eqcmdata)-1)[1:]:
        #while b <= len(eqcmdata)-2:
           # if eqcmdata['T'].loc[j] > df_cv['T'].iloc[b] == 1:
                #print 'test'
            if b >= len(eqcmdata)-2:  
                break
            else:
                b = b+1  
    eqcmdata_cv = df_cv
                
b = 1          
              
    
        #for k in cvheader[1:4]:
                # eqcmdata_cv = eqcmdata.iloc[a:b]
         #   eqcmdata_cv.loc[j,[k]] = df_cv.loc[b,[k]]
    #a = b        
             
            
            
    

#for i in range(len(cvlist)):
#    for j in range(5)[2:5]:
        
        
    
       



#df['T_new'] = map(zero, df)