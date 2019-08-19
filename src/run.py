#! /usr/bin/env python
#! -*- coding:utf-8 -*-
import time
import sys
import subprocess
import shutil
import os  

proc = subprocess.Popen(['python', 'app/python.py'],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
a = None
Test = None
while 1:
    text = proc.stdout.readline()[:-1] 
    if type(text) != str or text == '' and proc.poll() != None:
        break
    elif type(text) == str and len(text) > 6:
        Test = str(a)
        print text
if  Test == None:
    print "Check Error"
