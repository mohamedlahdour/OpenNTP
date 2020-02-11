#! /usr/bin/env python
#! -*- coding:utf-8 -*-
import time
import sys
import subprocess
import shutil
import os  
M00 = open('app/link/script00.py', "r" ).read() 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if M00 == 'CP1D':
    if os.path.exists('app/CP1D.so'):
        os.remove('app/CP1D.so')
    proc = subprocess.Popen(['f2py','-c','app/sources/TRANSPORT_1D_CP.f90','-m','CP1D'],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
    while 1:
        text = proc.stdout.readline()[:-1]
        if type(text) != str or text == '' and proc.poll() != None: 
            break
 
        elif type(text) == str and len(text) > 6:
            print text
    shutil.move('CP1D.so', 'app') 
    
elif M00 == 'SN1D':
    if os.path.exists('app/SlabSN.so'):
        os.remove('app/SlabSN.so') 
    proc = subprocess.Popen(['f2py','-c','app/sources/TRANSPORT_SN.f90','-m','SlabSN'], 
                                        stderr=subprocess.PIPE, 
                                        stdout=subprocess.PIPE)
    while 1:
        text = proc.stdout.readline()[:-1]
        if type(text) != str or text == '' and proc.poll() != None: 
            break
 
        elif type(text) == str and len(text) > 6:
            print text
    shutil.move('SlabSN.so', 'app')

elif M00 == 'SN2D':
    if os.path.exists('app/SlabSN2D.so'):
        os.remove('app/SlabSN2D.so') 

    proc = subprocess.Popen(['f2py','-c','app/sources/TRANSPORT_2D_SN.f90','-m','SlabSN2D'], 
                                        stderr=subprocess.PIPE, 
                                        stdout=subprocess.PIPE)
    while 1:
        text = proc.stdout.readline()[:-1]
        if type(text) != str or text == '' and proc.poll() != None: 
            break
 
        elif type(text) == str and len(text) > 6:
            print text
    shutil.move('SlabSN2D.so', 'app')

elif M00 == 'MOC1D':
    if os.path.exists('app/SlabMOC.so'):
        os.remove('app/SlabMOC.so')
    proc = subprocess.Popen(['f2py','-c','app/sources/TRANSPORT_MOC.f90','-m','SlabMOC'],
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
    while 1:
        text = proc.stdout.readline()[:-1]
        if type(text) != str or text == '' and proc.poll() != None: 
            break
 
        elif type(text) == str and len(text) > 6:
                print text
    shutil.move('SlabMOC.so', 'app')
else:
    print "select the calculation method."
   
