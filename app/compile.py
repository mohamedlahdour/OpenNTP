#! /usr/bin/python3
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
    myPopen = subprocess.Popen(['f2py3','-c','app/sources/TRANSPORT_CP.f90','-m','CP1D'], stdin = subprocess.PIPE,
                         stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = 'utf8')
    while 1:
        text = myPopen.stdout.readline()[:-1]
        if type(text) != str or text == '' and myPopen.poll() != None: 
            break
 
        elif type(text) == str and len(text) > 6:
            print (text)
    if os.path.exists('CP1D.cpython-38-x86_64-linux-gnu.so'):
        os.rename('CP1D.cpython-38-x86_64-linux-gnu.so', 'CP1D.so')
    shutil.move('CP1D.so', 'app')
    
elif M00 == 'SN1D':
    if os.path.exists('app/SlabSN.so'):
        os.remove('app/SlabSN.so') 
    myPopen = subprocess.Popen(['f2py3','-c','app/sources/TRANSPORT_SN.f90','-m','SlabSN'],
                                        stdin = subprocess.PIPE,
                                        stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = 'utf8')

    while 1:
        text = myPopen.stdout.readline()[:-1]
        if type(text) != str or text == '' and myPopen.poll() != None: 
            break
 
        elif type(text) == str and len(text) > 6:
            print (text)
    if os.path.exists('SlabSN.cpython-38-x86_64-linux-gnu.so'):
        os.rename('SlabSN.cpython-38-x86_64-linux-gnu.so', 'SlabSN.so')
    shutil.move('SlabSN.so', 'app')

elif M00 == 'SN2D':
    if os.path.exists('app/SN2D.so'):
        os.remove('app/SN2D.so')
    myPopen = subprocess.Popen(['f2py3','-c','app/sources/TRANSPORT_2D_SN.f90','-m','SN2D'],
                                        stderr=subprocess.PIPE, 
                                        stdout=subprocess.PIPE, encoding = 'utf8')

    while 1:
        text = myPopen.stdout.readline()[:-1]
        if type(text) != str or text == '' and myPopen.poll() != None: 
            break
 
        elif type(text) == str and len(text) > 6:
            print (text)
    if os.path.exists('SN2D.cpython-38-x86_64-linux-gnu.so'):
        os.rename('SN2D.cpython-38-x86_64-linux-gnu.so', 'SN2D.so')
    shutil.move('SN2D.so', 'app')

elif M00 == 'MOC1D':
    if os.path.exists('app/SlabMOC.so'):
        os.remove('app/SlabMOC.so')
    myPopen = subprocess.Popen(['f2py3','-c','app/sources/TRANSPORT_MOC.f90','-m','SlabMOC'],
                                        stdin = subprocess.PIPE,
                                        stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = 'utf8')
    while 1:
        text = myPopen.stdout.readline()[:-1]
        if type(text) != str or text == '' and myPopen.poll() != None: 
            break
 
        elif type(text) == str and len(text) > 6:
                print (text)
    if os.path.exists('SlabMOC.cpython-38-x86_64-linux-gnu.so'):
        os.rename('SlabMOC.cpython-38-x86_64-linux-gnu.so', 'SlabMOC.so')
    shutil.move('SlabMOC.so', 'app')
else:
    print ("select the calculation method.")
   
