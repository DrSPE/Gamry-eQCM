# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 10:53:57 2016

@author: se
"""

#Ordner: verschiebt einzelne Dateien in Ordner mit dem Namen der Datei -.DTA

import glob, os, shutil

folder = '/home/se/Dropbox/Masterarbeit/SEI/2016-08-11'

x=0

for file_path in glob.glob(os.path.join(folder, '*.*')):
    new_dir = file_path.rsplit('.', 1)[0]
    if file_path.rsplit('.', 1)[1] == 'DTA': 
        os.mkdir(os.path.join(folder, new_dir))
        shutil.move(file_path, os.path.join(new_dir, os.path.basename(file_path)))
        x = x+1
    
print ('directories created and files moved:')
print x