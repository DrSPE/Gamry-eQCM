# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 14:10:28 2016

@author: se
"""

#import of Gamry .DTA file with different paragraphs of CV data and a long QCM data paragraph in the end


#import glob, os
#os.chdir("/home/se/Dropbox/Masterarbeit/Projekt/data/raw/eQCM/C_beschichtet/#2/2016-05-19")
#for file in glob.glob("*.dta"):
f = open('/home/se/Dropbox/Masterarbeit/Projekt/data/raw/eQCM/C_beschichtet/#2/2016-05-19/QCM_RCVv8_10mV_all.DTA', 'r') #location of .DTA file
llist=[]

for line in f:
     line = line.strip()
     columns = line.split()
     llist.append(columns)
     llist = filter(None, llist) #removes empty lines
    
     cvlist=[]
     qcmlist=[]

   
for i in llist:
       if i[0].startswith('CURVE'):
                #print llist.index(i)
                cvlist.append(llist.index(i))
       if  'QCMDESCRIPTION' in i:
                 qcmlist = llist.index(i)
                 #print 'qcm position found:'           
                 #print qcmlist
                 #print llist.index(i)    
                        
                 #print 'finished'

print 'Es wurden'
print len(cvlist)
print 'CV Zyklen gefunden.'
print 'Diese starten in den Zeilen:'
print cvlist
print 'Die QCM Daten beginnen ab Zeile:'
print qcmlist


#Jetzt die Enzelnen Zyklen in einzelne Dateien schreiben? Mit for i in range(1, len(cvlist)+1)

