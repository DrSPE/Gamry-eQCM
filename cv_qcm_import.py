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
from matplotlib import ticker
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import griddata


"""
Enter data here:
e.g. active crystal area / Area between electrodes A
or if determined the shift of the QRE to Ferrocen

also check if the crystal constant Cq in the script was chosen correctly
for 10 MHz data it will be assumed that a 10MHz crystal from renlux was used
for 6MHz, that a 6MHz crystal from openQCM/novaetech it. was used
otherwise the constant will be calculated with the rounded resonant frequency, assuming a natural number in MHz
Cq =  (2* resonant_freq**2)/(math.sqrt(mu * p))  in [Hz cm²/ng] # here resonant_freq is in Hz =^= [MHz] *1e6
with            
mu = 2.947e20 #ng/(cm s^2)  ( Shear modulus of quartz for AT-cut crystal )
p = 2.648e9 #ng/cm^3        (density of quarz)

https://en.wikipedia.org/wiki/Sauerbrey_equation
"""
#Crystal data:
A = 0.205 #cm^2

#Ferrocen Kalibrierung in V:
U_k = 0.0

#Beschichtung:
carb = 'AC03'

#import of Gamry .DTA file with different paragraphs of CV data and a long QCM data paragraph in the end

#two options to enter the file path:

#1.
# directory path is split into folders: file(folder.dta) in folder in date in crystNr in direct
direct = '/home/se/Dropbox/Masterarbeit/'
date = '2016-08-11'
crystNr = 'SEI'
folder = 'eQCM-Elektrolyte2-CV_v1' 
filename = folder + '.DTA'
#data file is named as folder.dta !!!

u = direct +  crystNr +'/'+ date + '/'+ folder + '/' + filename #location of .DTA file

'''#2. option filepath is entered directly 
-> causes problems with the graphs since they refer to date!
u= '/home/se/Dropbox/Masterarbeit/Projekt/data/raw/eQCM/C_beschichtet/#2/2016-07-14/QCM_RCV_5mV_up_v1/QCM_RCV_5mV_up_v1.DTA'
'''
f = open(u, 'r') #open .DTA file
origin = os.path.dirname(u)+'/' #folder of .DTA file
# also possible: origin = direct + crystNr +'/'+ date + '/'+ folder + '/'
data=[] #initiate the list/tuple

for line in f:          #filles the list data with the data from .dta file
    line = line.strip()
    columns = line.split()
    data.append(columns)

data = filter(None, data) #removes empty lines in .dta file data
    
cvlist = [] #starting points of CV cycles/ row numbers in data list
qcmlist = [] #starting point of QCM data / row number in data list
qcmstart = [] #time the QCM was measuring before the electrochemical measurements were started and row postition of this information in data list

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

print ('Es wurden '+repr(len(cvlist)) +' CV Zyklen gefunden.')

if U_k ==0:     #in case no calibration was done, shows the voltage vs the reference used, in my case a Ag wire as Ag-QRE
    print 'Es wurde keine Kalibrierung der Referenz berücksichtigt'
    U_label = r' $ vs.$ $Ag-QRE $ '
else:           #in case a calibration was done (shift has to be entered at the top of this script!!!) recalculates the voltage data with this shift and shows voltage vs half wave potential of the Ferrocen redox couple
    print 'Es wurde eine Verschiebung der Ag-QRE von {0} V berücksichtigt und E(1/2)[Ferrocen] = 0 V gesetzt '.format(U_k)
    U_label = r' $vs.$ $E_{1/2}(Ferrocen)$ '
    
    
header = data[0:(cvlist[0])] #includes the text in the header of the .dta file

#CV Zyklen splitten:

cvname = [] #initiate cv data names

a=0 
for i in cvlist:
    a = a+1
    cvname.append('CV_{0:02d}'.format(a))
    
cvdata = []    #initiate echem data collection, enthält alle CV Zyklen (getrennt)

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
          
#takes the value in the first row to calculate delta f_s and delta f_p
fs0 = qcmdata[2][2]
fp0 = qcmdata[2][3]

#saving QCM and CV data to files


#save header

nameheader = origin + 'header.txt'
file = open(nameheader, 'w')
wr = csv.writer(file, quoting=csv.QUOTE_ALL)
wr.writerows((header))
file.close()


#saves each CV cycle data to a seperate file
for i in range(len(cvname)):
    namecv = origin + cvname[i] + '.txt'
    cv_loop =np.asarray(cvdata[i][2:])
    cv_loop = np.asarray(cv_loop[:][:,1:5], dtype = np.float32)
    np.savetxt(namecv, cv_loop, delimiter=',', header='T/s, Q/C, U_ref/V, I/A')
    

cvtime = [] #starting time of each cv cycle
for i in range(len(cvname)):
    cvtime.append(float(cvdata[i][2][1]))



#modify QCM data to start at 0 seconds    
start = 0 

while start < len(qcmdata):
      if start <2:          #skips the first 2 row, since they contain the heding and unit
          start= start + 1
      else:
          qcmdata[start][1] = float(qcmdata[start][1]) - qcm_time_zero
          start = start + 1
     
#saves the qcmdata with modified time to file (all in one)  
nameqcm = origin + 'qcm_data.txt'
file = open(nameqcm, 'w')
for i in range(len(qcmdata)):
  file.write("%s \n" % qcmdata[i])
file.close()

# End of wrtining original data to seperate files

#creates file names to be used for the seperated qcmdata later on
q = 0
qcmname = []
for i in cvname:
   q = q+1
   qcmname.append('QCM_{0:02d}'.format(q))

#saves the qcmdata as a file (all in one) without the unnecessary parts
eqcmdata = np.asarray(qcmdata)[2:][:,1:4]
eqcmdata = np.asarray(eqcmdata, dtype=np.float32)
eqcmname = origin + 'eqcm_data.txt'
np.savetxt(eqcmname, eqcmdata, delimiter=',', header='T/s , fs/MHz , fp/MHz')

#split/match qcmdata to CV cycles and saves the location of each new cycle to the list qcmtime
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

#saves the qcmdata matching each CV cacle to seperate files
l = 0
for i in range(len(qcmtime)):
    if qcmtime[i]<qcmtime[-1]:
        np.savetxt(origin + qcmname[i] + '.txt', eqcmdata[l:qcmtime[i+1]], delimiter=',', header='T/s , fs/MHz , fp/MHz')
        l = qcmtime[i+1]
    else:
        np.savetxt(origin +qcmname[i] + '.txt', eqcmdata[qcmtime[-1]:], delimiter=',', header='T/s , fs/MHz , fp/MHz')
        
'''
imports the data from one cycle per loop, then interpolates the CV data since they have 
a bigger distance between each data point, with interp1d (linear interpolation of scipy):fit_...
then they make sure that the data for the interpolation (QCM =^=time data of qcm measurement)
starts after the fitted area (1.while loop) and ends beofre the fitted area (2. loop), otherwise it would be extrapolation and give an error
then the interpolation is done: ..._new
'''
for i in range(len(cvname)):
    time = np.loadtxt(origin + cvname[i] + '.txt', delimiter=',', skiprows=1, usecols =(0,))
    Q = np.loadtxt(origin + cvname[i] + '.txt', delimiter=',', skiprows=1, usecols =(1,))
    U = np.loadtxt(origin + cvname[i] + '.txt', delimiter=',', skiprows=1, usecols =(2,))
    I = np.loadtxt(origin + cvname[i] + '.txt', delimiter=',', skiprows=1, usecols =(3,))
    QCM = np.loadtxt(origin + qcmname[i] + '.txt', delimiter=',', skiprows=1, usecols =(0,))
    fs = np.loadtxt(origin + qcmname[i] + '.txt', delimiter=',', skiprows=1, usecols =(1,))
    fp = np.loadtxt(origin + qcmname[i] + '.txt', delimiter=',', skiprows=1, usecols =(2,))
    fitQ = interp1d(time, Q)
    fitU = interp1d(time, U)
    fitI = interp1d(time, I)
    while time[0]>QCM[0]:
        QCM = QCM[1:]
        fs  =  fs[1:]
        fp  =  fp[1:]
    while time[-1]<QCM[-1]:
        QCM = QCM[:-1]
        fs  =  fs[:-1]
        fp  =  fp[:-1]
    Q_new = fitQ(QCM)
    Q_new = Q_new * 1E3
    U_new = fitU(QCM)
    U_new = U_new - U_k
    I_new = fitI(QCM)
    I_new = I_new * 1E6 #the current is converted from A to nA 
    dfs = (fs- float(fs0))*1E6  #calculating delta f_s
    dfp = (fp- float(fp0))*1E6  #and delta f_p
    dfsp = fp*1E6 - fs*1E6
    
    '''
    chooses the crystal constant in the first loop
    '''
    if i ==0:   
        if int(round(np.mean(fs), 0)) == 6:
            Cq = 0.0812 #Hz cm²/ng
            print 'Die Massenänderung wurde für einen 6 MHz Kristall von novaetech berechnet.'
        if int(round(np.mean(fs), 0)) == 10:
            Cq = 0.223  #Hz cm²/ng
            print 'Die Massenänderung wurde für einen 10 MHz Kristall von renlux berechnet.'
        if int(round(np.mean(fs), 0)) != 6:
            if int(round(np.mean(fs), 0)) != 10:
                resonant_freq = int(round(np.mean(fs), 0)) * 1e6
                mu = 2.947e20 #ng/(cm s^2)
                p = 2.648e9 #ng/cm^3
            
                Cq =  (2* resonant_freq**2)/(math.sqrt(mu * p)) #Hz cm²/ng'''
                print 'Die Massenänderung wurde für einen nicht bekannten {} MHz Kristall berechnet. Unbedingt prüfen ob die Parameter wie Fläche etc. stimmen!'.format(int(round(np.mean(fs), 0)))
        print 'Ist die Fläche von A = {} cm² korrekt?'.format(A)
    
        '''
        caclulates the mass change/area (dm1) and the absolute mass change (dm2) 
        and some kind of quality factor wf which should not change greatly during the measurement 
        otherwise there might be a viscoelastic behaviour
        '''    
    
    dm1 = -fs/Cq    
    dm2 = dm1*A
    
    wf = (fs-fp)/((fs+fp)/2)
    
    '''
    saves the eQCM data to seperate files for each cycle and then plots a few figures and saves them
    '''    
    np.savetxt(origin + 'eQCM_{0:02d}.txt'.format(i+1), np.c_[QCM, fs, fp, dfs, dfp,dm1, dm2, Q_new, U_new, I_new], delimiter=',', newline='\n', header='T/s , fs/MHz , fp/MHz, dfs, dfp, dm1/[ng/cm²], dm2/ng, Q/mC, U_Ref/V, I/µA')
    
    
    
   
    fig1 = plt.figure()
    fig1.suptitle(date + ': ' + folder + ' cycle {0:02d}'.format(i+1) )
    
     
    ax1 = fig1.add_subplot(221)
    ax1.plot(Q_new, dm2, 'g-')
    
    ax1.set_xlabel(r'$Q / \mathrm{mC}$')
    ax1.set_ylabel(r'$\Delta m / \mathrm{ng}$')
    #ax1.set_title('Crystal ' + crystNr + ' coated with ' + carb + ', measured on ' + date + ': cycle {0:02d}'.format(i+1) )
    
    
    ax2 = fig1.add_subplot(222)
    ax2.plot(QCM, dfs, 'b-')
    ax2.set_ylabel(r'$ \Delta f_s / \mathrm{Hz}$') 
    ax2.set_xlabel(r'$t / \mathrm{s}$')
    
    ax3 = fig1.add_subplot(223)
    ax3.plot(QCM, dfsp, 'g-')
    ax3.set_ylabel(r'$f_p - f_s / \mathrm{Hz}$')   
    #ax3.set_ylabel(r'$\Delta W/ \Delta f$')
    ax3.set_xlabel(r'$t / \mathrm{s}$')
    
    ax4 = fig1.add_subplot(224)
    ax4.plot(U_new, I_new, 'g-')
    ax4.set_ylabel(r'$I / \mathrm{\mu A}$')
    ax4.set_xlabel(r'$U$' + U_label + r'$/ \mathrm{V}$')
    
    fig1.tight_layout(pad=2.0)  

    
    fig1.savefig(origin + 'cycle_{0:02d}.pdf'.format(i+1), bbox_inches='tight')
    #plt.clf()
    plt.close()   
    
  

    
    
    '''
    fig = plt.figure(1)
    plt.subplot(221)
    plt.plot(Q_new, dm2, 'g-')
    plt.ylabel(r'$\Delta m / \mathrm{ng}$')
    plt.xlabel(r'$Q / \mathrm{mC}$')
    
    plt.subplot(222)
    plt.plot(U_new, dm2, 'g-')
    plt.ylabel(r'$\Delta m / \mathrm{ng}$')
    plt.xlabel(r'$U$' + U_label + r'$/ \mathrm{V}$')
    
    plt.subplot(223)
    plt.plot(QCM, wf, 'g-')
    plt.ylabel(r'$\Delta W/ \Delta f$')
    plt.xlabel(r'$t / \mathrm{s}$')
    
    plt.subplot(224)
    plt.plot(U_new, I_new, 'g-')
    plt.ylabel(r'$I / \mathrm{\mu A}$')
    plt.xlabel(r'$U$' + U_label + r'$/ \mathrm{V}$')
    
    fig.tight_layout()

    
    savefig(origin + 'cycle_{0:02d}.pdf'.format(i+1), bbox_inches='tight')
    #plt.clf()
    plt.close()   
    '''
    
    
    
    """
    Vielleicht noch eine Kurve mit Ladung in e- als x-Achse 
    und u(elementary mass unit 12C/12) als y-Achse?
    
    Für die Bilder Crystal und Datum als Überschrift, und Beschichtung?
    
    """
   
   
   
   
   
   
   
   
   
   
   