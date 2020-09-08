#! /usr/bin/ python3
#! -*- coding:utf-8 -*-
import sys
import os
from multiprocessing import Queue
from PyQt5 import QtCore, QtGui, QtWidgets
from app.mainwindow import Ui_MainWindow
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
import subprocess
import sys
import webbrowser
import numpy as np
import numpy 
from subprocess import Popen, STDOUT, PIPE
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import json
import math 
import matplotlib.patches as mpatch
from decimal import Decimal
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse, Arc
from matplotlib.patches import Circle, PathPatch
from matplotlib.transforms import Bbox
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from functools import partial
from matplotlib.widgets import Button
import matplotlib.pyplot
import time
# The new Stream Object which replaces the default stream associated with sys.stdout
# This object just puts data in a queue!
class EmittingStream(QtCore.QObject):
    textWritten = QtCore.pyqtSignal(str)
 
    def write(self,text):
        self.textWritten.emit(str(text))
        pass

    def flush(self):
        pass


# An example QObject (to be run in a QThread) which outputs information with print
class LongRunningThing(QObject):
    @pyqtSlot()
    def run(self):
        pass
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window1(QWidget):
    def __init__(self,npc,nxpc,parent=None): 
        super(Window1, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)

        self.nx=nxpc
        self.lab4  = [0]*nxpc
        self.lineEdit8 = [0]*nxpc*npc

        self.lab6  = [0]*nxpc
        self.lineEdit10 = [0]*nxpc*npc
        self.lab7  = [0]*nxpc
        self.lineEdit11 = [0]*nxpc*npc
        self.dx  = [0]*nxpc*npc

        self.Vect6 = [0]*nxpc*npc
        self.Vect7 = [0]*nxpc*npc
        self.regmat = [[[]]]
        self.nfmesh_x = [[[]]]
        self.types =  []
        M01 = open('app/link/script01.py', "r" ).read() 
        for i in range(npc):
            self.types.append("Pin Cell %s" %(i+1))
        vbox = QVBoxLayout()
        tabWidget = QTabWidget()
        num=0
        for name in self.types:
            TabContact = QWidget()
            tabWidget.setFont(QtGui.QFont("Sanserif", 10))
            tabWidget.addTab(TabContact, name[0:11])
            vbox.addWidget(tabWidget)
            self.setLayout(vbox)


            if M01 == 'Slab Geometry':
                groupBox = QGroupBox("Size for each region per [cm]")
            else:
                groupBox = QGroupBox("Ray for each annular geometry per [cm]")
            
            groupBox.setStyleSheet('QGroupBox:title {color: blue;}')
            typetablayout = QGridLayout()
            for i in range(nxpc):
                if M01 == 'Slab Geometry':
                    self.lab4[i] = QLabel("X %s" %(i+1))
                else:
                    self.lab4[i] = QLabel("R %s" %(i+1))
                self.lab4[i].setAlignment(Qt.AlignCenter)
                typetablayout.addWidget(self.lab4[i], 1, i+1)
                self.lineEdit8[nxpc*num+i] = QLineEdit()
                typetablayout.addWidget(self.lineEdit8[nxpc*num+i], 2, i+1)
                self.lineEdit8[nxpc*num+i].insert(str(self.dx[nxpc*num+i]))

            groupBox.setLayout(typetablayout)  

            if M01 == 'Slab Geometry':
                groupBox2 = QGroupBox("Which material fills each region")
            else:
                groupBox2 = QGroupBox("Which material fills each annular geometry")
            groupBox2.setStyleSheet('QGroupBox:title {color: blue;}')
            typetablayout2 = QGridLayout()

            for i in range(nxpc):
                if M01 == 'Slab Geometry':
                    self.lab6[i] = QLabel("X %s" %(i+1))
                else:
                    self.lab6[i] = QLabel("R %s" %(i+1))
                self.lab6[i].setAlignment(Qt.AlignCenter)
                typetablayout2.addWidget(self.lab6[i], 6, i+1)
                self.lineEdit10[nxpc*num+i] = QLineEdit()
                typetablayout2.addWidget(self.lineEdit10[nxpc*num+i], 7, i+1)
                self.lineEdit10[nxpc*num+i].insert(str(self.Vect6[nxpc*num+i]))
            groupBox2.setLayout(typetablayout2)


            if M01 == 'Slab Geometry':
                groupBox3 = QGroupBox("Number of fine meshes per region")
            else:
                groupBox3 = QGroupBox("Number of fine meshes per annular geometry")
            groupBox3.setStyleSheet('QGroupBox:title {color: blue;}')
            typetablayout3 = QGridLayout()

            for i in range(nxpc):
                if M01 == 'Slab Geometry':
                    self.lab7[i] = QLabel("X %s" %(i+1))
                else:
                    self.lab7[i] = QLabel("R %s" %(i+1))
                self.lab7[i].setAlignment(Qt.AlignCenter)
                typetablayout3.addWidget(self.lab7[i], 9, i+1)
                self.lineEdit11[nxpc*num+i] = QLineEdit()
                typetablayout3.addWidget(self.lineEdit11[nxpc*num+i], 10, i+1)
                self.lineEdit11[nxpc*num+i].insert(str(self.Vect7[nxpc*num+i]))

            groupBox3.setLayout(typetablayout3)
            mainLayout = QVBoxLayout(TabContact)
            mainLayout.addWidget(groupBox)
            mainLayout.addWidget(groupBox2)
            mainLayout.addWidget(groupBox3)

            #réer un Bouton
            if num == (len(self.types)-1):
                self.bouton = QPushButton(u"Save and Colse")
                self.bouton.clicked.connect(partial(self.save6,num))
                vbox.addWidget(self.bouton)

            num+=1
    def save6(self,num):
	#connexion avec la fenetre main
        del self.dx[:]
        for k in range(len(self.types)):
            self.dx.append([])
            for i in range(self.nx):
                self.dx[k].append(eval(self.lineEdit8[self.nx*k+i].text()))   

        #------------------------------
        del self.regmat[:]
        for k in range(len(self.types)):
            m=0
            self.regmat.append([])
            for i in range(self.nx):
                self.regmat[k].append(eval(self.lineEdit10[self.nx*k+m].text()))
                m+=1
        #------------------------------
        del self.nfmesh_x[:]
        for k in range(len(self.types)):
            m=0
            self.nfmesh_x.append([])
            for i in range(self.nx):
                self.nfmesh_x[k].append(eval(self.lineEdit11[self.nx*k+m].text()))
                m+=1

        if num == (len(self.types)-1):
            self.close() 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window2(QWidget):
    def __init__(self,npca,na,parent=None):
        super(Window2, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        self.setMinimumSize(QtCore.QSize(100, 100))
        self.npca = npca
        self.lab6  = [0]*npca*na
        self.lineEdit10 = [0]*npca*na
        self.Vect6 = [0]*npca*na
        self.assembly = []
        self.layout = QVBoxLayout(self)
        typetab = QTabWidget(self)  
        typetab.setFont(QtGui.QFont("Sanserif", 10)) 
        self.types =  []
        for i in range(na):
            self.types.append("Assembly %s" %(i+1))
        num=0
        for name in self.types:
            tab =  QWidget()
            typetab.addTab(tab, name[0:10])
            typetablayout = QGridLayout(tab)

            for i in range(npca):
                self.lab6[i] = QLabel("%s" %(i+1))
                self.lab6[i].setAlignment(Qt.AlignCenter)
                typetablayout.addWidget(self.lab6[i], 5, i+1)
                self.lineEdit10[npca*num+i] = QLineEdit()
                typetablayout.addWidget(self.lineEdit10[npca*num+i], 6, i+1)
                self.lineEdit10[npca*num+i].insert(str(self.Vect6[npca*num+i]))

            self.layout.addWidget(typetab)
            self.setLayout(self.layout)
            #réer un Bouton
            if num == (len(self.types)-1):
                self.bouton = QPushButton(u"Save and Colse")
                self.bouton.clicked.connect(partial(self.save7,num))
                self.layout.addWidget(self.bouton)
            num+=1
    def save7(self,num):
	#connexion avec la fenetre main
        del self.assembly[:]
        for k in range(len(self.types)):
            m=0
            self.assembly.append([])
            for i in range(self.npca):
                self.assembly[k].append(eval(self.lineEdit10[self.npca*k+m].text()))
                m+=1
        if num == (len(self.types)-1):
            self.close() 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window3(QWidget):
    def __init__(self,na,parent=None):
        super(Window3, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        self.setMinimumSize(QtCore.QSize(100, 100))
        self.na = na
        self.lab4  = [0]*na
        self.lineEdit8 = [0]*na
        self.core   = [0]*na
        self.layout = QVBoxLayout(self)    
        typetablayout = QGridLayout()
        #Main group box
        box = QGroupBox()
        #self.main_group_box.setStyleSheet("QGroupBox{font-size: 10px}")
        box.setTitle("Core geometry:")
        box.setStyleSheet('QGroupBox:title {color: blue;}')
        #self.main_group_box.setLayout(self.layout)

        for i in range(na):
            self.lab4[i] = QLabel("%s" %(i+1))
            self.lab4[i].setAlignment(Qt.AlignCenter)
            typetablayout.addWidget(self.lab4[i], 6, i+1)
            self.lineEdit8[i] = QLineEdit()
            typetablayout.addWidget(self.lineEdit8[i], 7, i+1)
            self.lineEdit8[i].insert(str(self.core[i]))


        box.setLayout(typetablayout)
        self.layout.addWidget(box)
        self.setLayout(self.layout)

        #réer un Bouton
        self.bouton = QPushButton(u"Save and Colse")
        self.bouton.clicked.connect(self.save3)
        self.layout.addWidget(self.bouton)


    def save3(self):
	#connexion avec la fenetre main
        del self.core[:] 
        for i in range(self.na):
            self.core.append(eval(self.lineEdit8[i].text()))   
        self.close()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window9(QWidget):
    def __init__(self,ngroup,nmat,parent=None):
        super(Window9, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        self.setMinimumSize(QtCore.QSize(100, 100))

        self.ng = ngroup
        self.nm = nmat
        self.SigT  = [0]*ngroup*nmat
        self.tab2  = [[0]*ngroup]*nmat
        self.lineEdit4 = [0]*ngroup*nmat 
        lab4  = [0]*ngroup*nmat

        self.layout = QVBoxLayout(self)    
        typetablayout = QGridLayout()
        #Main group box
        box = QGroupBox()
        #self.main_group_box.setStyleSheet("QGroupBox{font-size: 10px}")
        box.setTitle("Total Cross Section (SigT)")
        box.setStyleSheet('QGroupBox:title {color: blue;}')
        #self.main_group_box.setLayout(self.layout)

        m = 0
        for j in range(nmat):
            lab4[j] = QLabel("Material %s" %(j+1))
            lab4[j].setAlignment(Qt.AlignCenter)
            typetablayout.addWidget(lab4[j], j+2, 0)
            for i in range(ngroup):
                lab4[i] = QLabel("G %s" %(i+1))
                lab4[i].setAlignment(Qt.AlignCenter)
                typetablayout.addWidget(lab4[i], 1, i+1)
                self.lineEdit4[m] = QLineEdit()
                typetablayout.addWidget(self.lineEdit4[m], j+2, i+1)
                self.lineEdit4[m].insert(str(self.SigT[m]))
                m=m+1

        box.setLayout(typetablayout)
        self.layout.addWidget(box)
        self.setLayout(self.layout)

        #réer un Bouton
        self.bouton = QPushButton(u"Save and Colse")
        self.bouton.clicked.connect(self.save8)
        self.layout.addWidget(self.bouton)

    def save8(self):
        m=0
        self.SigT = []
        for i in range(self.nm):
            self.SigT.append([])
            for j in range(self.ng):
                self.SigT[i].append(eval(self.lineEdit4[m].text()))
                m+=1
        self.close()  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window10(QWidget):
    def __init__(self,ngroup,nmat,parent=None):
        super(Window10, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        self.setMinimumSize(QtCore.QSize(100, 100))

        self.ng = ngroup
        self.nm = nmat
        self.SigF  = [0]*ngroup*nmat
        self.tab2  = [[0]*ngroup]*nmat
        self.lineEdit4 = [0]*ngroup*nmat 
        lab4  = [0]*ngroup*nmat

        self.layout = QVBoxLayout(self)    
        typetablayout = QGridLayout()
        #Main group box
        box = QGroupBox()
        #self.main_group_box.setStyleSheet("QGroupBox{font-size: 10px}")
        box.setTitle("Fission Cross Section (SigF)")
        box.setStyleSheet('QGroupBox:title {color: blue;}')
        #self.main_group_box.setLayout(self.layout)

        m = 0
        for j in range(nmat):
            lab4[j] = QLabel("Material %s" %(j+1))
            lab4[j].setAlignment(Qt.AlignCenter)
            typetablayout.addWidget(lab4[j], j+2, 0)
            for i in range(ngroup):
                lab4[i] = QLabel("G %s" %(i+1))
                lab4[i].setAlignment(Qt.AlignCenter)
                typetablayout.addWidget(lab4[i], 1, i+1)
                self.lineEdit4[m] = QLineEdit()
                typetablayout.addWidget(self.lineEdit4[m], j+2, i+1)
                self.lineEdit4[m].insert(str(self.SigF[m]))
                m=m+1

        box.setLayout(typetablayout)
        self.layout.addWidget(box)
        self.setLayout(self.layout)

        #réer un Bouton
        self.bouton = QPushButton(u"Save and Colse")
        self.bouton.clicked.connect(self.save8)
        self.layout.addWidget(self.bouton)

    def save8(self):
        m=0
        self.SigF = []
        for i in range(self.nm):
            self.SigF.append([])
            for j in range(self.ng):
                self.SigF[i].append(eval(self.lineEdit4[m].text()))
                m+=1
        self.close()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window11(QWidget):
    def __init__(self,ngroup,nmat,parent=None):
        super(Window11, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        self.setMinimumSize(QtCore.QSize(100, 100))

        self.ng = ngroup
        self.nm = nmat
        self.NuSigF  = [0]*ngroup*nmat
        self.tab2  = [[0]*ngroup]*nmat
        self.lineEdit4 = [0]*ngroup*nmat 
        lab4  = [0]*ngroup*nmat

        self.layout = QVBoxLayout(self)    
        typetablayout = QGridLayout()
        #Main group box
        box = QGroupBox()
        #self.main_group_box.setStyleSheet("QGroupBox{font-size: 10px}")
        box.setTitle("NuFission Cross Section (NuSigF)")
        box.setStyleSheet('QGroupBox:title {color: blue;}')
        #self.main_group_box.setLayout(self.layout)

        m = 0
        for j in range(nmat):
            lab4[j] = QLabel("Material %s" %(j+1))
            lab4[j].setAlignment(Qt.AlignCenter)
            typetablayout.addWidget(lab4[j], j+2, 0)
            for i in range(ngroup):
                lab4[i] = QLabel("G %s" %(i+1))
                lab4[i].setAlignment(Qt.AlignCenter)
                typetablayout.addWidget(lab4[i], 1, i+1)
                self.lineEdit4[m] = QLineEdit()
                typetablayout.addWidget(self.lineEdit4[m], j+2, i+1)
                self.lineEdit4[m].insert(str(self.NuSigF[m]))
                m=m+1

        box.setLayout(typetablayout)
        self.layout.addWidget(box)
        self.setLayout(self.layout)

        #réer un Bouton
        self.bouton = QPushButton(u"Save and Colse")
        self.bouton.clicked.connect(self.save8)
        self.layout.addWidget(self.bouton)

    def save8(self):
        m=0
        self.NuSigF = []
        for i in range(self.nm):
            self.NuSigF.append([])
            for j in range(self.ng):
                self.NuSigF[i].append(eval(self.lineEdit4[m].text()))
                m+=1
        self.close() 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window12(QWidget):
    def __init__(self,ngroup,nmat,parent=None):
        super(Window12, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        self.setMinimumSize(QtCore.QSize(100, 100))
        self.ng = ngroup
        self.nm = nmat
        #self.ord = order
        self.lab6  = [0]*ngroup*ngroup*nmat
        self.lineEdit10 = [0]*ngroup*ngroup*nmat
        self.SigS = [0]*ngroup*ngroup*nmat
        self.layout = QVBoxLayout(self)
        typetab = QTabWidget(self)  
        typetab.setFont(QtGui.QFont("Sanserif", 10)) 
        self.types =  []

        self.lab6[0] = QLabel("<font color=blue > Scatter Matrix Cross Section (SigS)</font>")
        self.layout.addWidget(self.lab6[0])

        for i in range(nmat):
            self.types.append("Material %s" %(i+1))
        num=0
        for name in self.types:
            tab =  QWidget()
            typetab.addTab(tab, name[0:10])
            typetablayout = QGridLayout(tab)
            m = 0
            for j in range(ngroup):
                self.lab6[j] = QLabel("G %s" %(j+1))
                self.lab6[j].setAlignment(Qt.AlignCenter)
                typetablayout.addWidget(self.lab6[j], j+7, 0)
                for i in range(ngroup):
                    self.lab6[i] = QLabel("G %s" %(i+1))
                    self.lab6[i].setAlignment(Qt.AlignCenter)
                    typetablayout.addWidget(self.lab6[i], 6, i+1)
                    self.lineEdit10[ngroup*ngroup*num+m] = QLineEdit()
                    typetablayout.addWidget(self.lineEdit10[ngroup*ngroup*num+m], j+7, i+1)
                    self.lineEdit10[ngroup*ngroup*num+m].insert(str(self.SigS[ngroup*ngroup*num+m]))
                    m+=1

            self.layout.addWidget(typetab)
            self.setLayout(self.layout)
            #réer un Bouton
            if num == (len(self.types)-1):
                self.bouton = QPushButton(u"Save and Colse")
                self.bouton.clicked.connect(partial(self.save7,num))
                self.layout.addWidget(self.bouton)
            num+=1
    def save7(self,num):
	#connexion avec la fenetre main
        del self.SigS[:]
        for k in range(len(self.types)):
            m=0
            self.SigS.append([])
            for j in range(self.ng):
                self.SigS[k].append([])
                for i in range(self.ng):
                    self.SigS[k][j].append(eval(self.lineEdit10[self.ng*self.ng*k+m].text()))
                    m+=1
        if num == (len(self.types)-1):
            self.close() 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window13(QWidget):
    def __init__(self,ngroup,nmat,parent=None):
        super(Window13, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        self.setMinimumSize(QtCore.QSize(100, 100))

        self.ng = ngroup
        self.nm = nmat
        self.Chi  = [0]*ngroup*nmat
        self.tab2  = [[0]*ngroup]*nmat
        self.lineEdit4 = [0]*ngroup*nmat 
        lab4  = [0]*ngroup*nmat

        self.layout = QVBoxLayout(self)    
        typetablayout = QGridLayout()
        #Main group box
        box = QGroupBox()
        #self.main_group_box.setStyleSheet("QGroupBox{font-size: 10px}")
        box.setTitle("Density Function for Neutrons (Chi)")
        box.setStyleSheet('QGroupBox:title {color: blue;}')
        #self.main_group_box.setLayout(self.layout)

        m = 0
        for j in range(nmat):
            lab4[j] = QLabel("Material %s" %(j+1))
            lab4[j].setAlignment(Qt.AlignCenter)
            typetablayout.addWidget(lab4[j], j+2, 0)
            for i in range(ngroup):
                lab4[i] = QLabel("G %s" %(i+1))
                lab4[i].setAlignment(Qt.AlignCenter)
                typetablayout.addWidget(lab4[i], 1, i+1)
                self.lineEdit4[m] = QLineEdit()
                typetablayout.addWidget(self.lineEdit4[m], j+2, i+1)
                self.lineEdit4[m].insert(str(self.Chi[m]))
                m=m+1

        box.setLayout(typetablayout)
        self.layout.addWidget(box)
        self.setLayout(self.layout)

        #réer un Bouton
        self.bouton = QPushButton(u"Save and Colse")
        self.bouton.clicked.connect(self.save8)
        self.layout.addWidget(self.bouton)

    def save8(self):
        m=0
        self.Chi = []
        for i in range(self.nm):
            self.Chi.append([])
            for j in range(self.ng):
                self.Chi[i].append(eval(self.lineEdit4[m].text()))
                m+=1
        self.close()  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window6(QWidget):
    def __init__(self,nx,ny,np,parent=None):
        super(Window6, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        #self.setMinimumSize(QtCore.QSize(100, 400))
        self.nx=nx
        self.ny=ny
        self.lab4  = [0]*nx
        self.lineEdit8 = [0]*nx*np
        self.lab5  = [0]*nx
        self.lineEdit9 = [0]*ny*np
        self.lab6  = [0]*nx*ny
        self.lineEdit10 = [0]*nx*ny*np
        self.lab7  = [0]*nx*ny
        self.lineEdit11 = [0]*nx*ny*np
        self.dx  = [0]*nx*np
        self.dy  = [0]*ny*np
        self.Vect6 = [0]*nx*ny*np
        self.Vect7 = [0]*nx*ny*np
        self.regmat = [[[]]]
        self.nfmesh_xy = [[[]]]
        self.types =  []
        for i in range(np):
            self.types.append("Pin Cell %s" %(i+1))
        vbox = QVBoxLayout()
        tabWidget = QTabWidget()
        num=0
        for name in self.types:
            TabContact = QWidget()
            tabWidget.setFont(QtGui.QFont("Sanserif", 10))
            tabWidget.addTab(TabContact, name[0:11])
            vbox.addWidget(tabWidget)
            self.setLayout(vbox)

            groupBox = QGroupBox("Size for each region per [cm]")
            groupBox.setStyleSheet('QGroupBox:title {color: blue;}')
            typetablayout = QGridLayout()
            for i in range(nx):
                self.lab4[i] = QLabel("X %s" %(i+1))
                self.lab4[i].setAlignment(Qt.AlignCenter)
                typetablayout.addWidget(self.lab4[i], 1, i+1)
                self.lineEdit8[nx*num+i] = QLineEdit()
                typetablayout.addWidget(self.lineEdit8[nx*num+i], 2, i+1)
                self.lineEdit8[nx*num+i].insert(str(self.dx[nx*num+i]))
            for j in range(ny):
                self.lab5[j] = QLabel("Y %s" %(j+1))
                self.lab5[j].setAlignment(Qt.AlignCenter)
                typetablayout.addWidget(self.lab5[j], 3, j+1)
                self.lineEdit9[ny*num+j] = QLineEdit()
                typetablayout.addWidget(self.lineEdit9[ny*num+j], 4, j+1)
                self.lineEdit9[ny*num+j].insert(str(self.dy[ny*num+j]))
            groupBox.setLayout(typetablayout)

            groupBox2 = QGroupBox("Which material fills each region")
            groupBox2.setStyleSheet('QGroupBox:title {color: blue;}')
            typetablayout2 = QGridLayout()
            m = 0
            for j in range(ny):
                self.lab6[j] = QLabel("Y %s" %(j+1))
                self.lab6[j].setAlignment(Qt.AlignCenter)
                typetablayout2.addWidget(self.lab6[j], j+7, 0)
                for i in range(nx):
                    self.lab6[i] = QLabel("X %s" %(i+1))
                    self.lab6[i].setAlignment(Qt.AlignCenter)
                    typetablayout2.addWidget(self.lab6[i], 6, i+1)
                    self.lineEdit10[nx*ny*num+m] = QLineEdit()
                    typetablayout2.addWidget(self.lineEdit10[nx*ny*num+m], j+7, i+1)
                    self.lineEdit10[nx*ny*num+m].insert(str(self.Vect6[nx*ny*num+m]))
                    m+=1
            groupBox2.setLayout(typetablayout2)

            groupBox3 = QGroupBox("Number of fine meshes per region")
            groupBox3.setStyleSheet('QGroupBox:title {color: blue;}')
            typetablayout3 = QGridLayout()
            m = 0
            for j in range(ny):
                self.lab7[j] = QLabel("Y %s" %(j+1))
                self.lab7[j].setAlignment(Qt.AlignCenter)
                typetablayout3.addWidget(self.lab7[j], 10+ny+j, 0)
                for i in range(nx):
                    self.lab7[i] = QLabel("X %s" %(i+1))
                    self.lab7[i].setAlignment(Qt.AlignCenter)
                    typetablayout3.addWidget(self.lab7[i], 9+ny, i+1)
                    self.lineEdit11[nx*ny*num+m] = QLineEdit()
                    typetablayout3.addWidget(self.lineEdit11[nx*ny*num+m], j+10+ny, i+1)
                    self.lineEdit11[nx*ny*num+m].insert(str(self.Vect7[nx*ny*num+m]))
                    m=m+1
            groupBox3.setLayout(typetablayout3)
            mainLayout = QVBoxLayout(TabContact)
            mainLayout.addWidget(groupBox)
            mainLayout.addWidget(groupBox2)
            mainLayout.addWidget(groupBox3)

            #réer un Bouton
            if num == (len(self.types)-1):
                self.bouton = QPushButton(u"Save and Colse")
                self.bouton.clicked.connect(partial(self.save6,num))
                vbox.addWidget(self.bouton)

            num+=1
    def save6(self,num):
	#connexion avec la fenetre main
        del self.dx[:]
        del self.dy[:] 
        for k in range(len(self.types)):
            self.dx.append([])
            for i in range(self.nx):
                self.dx[k].append(eval(self.lineEdit8[self.nx*k+i].text()))   
            self.dy.append([])      
            for j in range(self.ny):
                self.dy[k].append(eval(self.lineEdit9[self.ny*k+j].text()))
        #------------------------------
        del self.regmat[:]
        for k in range(len(self.types)):
            m=0
            self.regmat.append([])
            for j in range(self.ny):
                self.regmat[k].append([])
                for i in range(self.nx):
                    self.regmat[k][j].append(eval(self.lineEdit10[self.nx*self.ny*k+m].text()))
                    m+=1
        #------------------------------
        del self.nfmesh_xy[:]
        for k in range(len(self.types)):
            m=0
            self.nfmesh_xy.append([])
            for j in range(self.ny):
                self.nfmesh_xy[k].append([])
                for i in range(self.nx):
                    self.nfmesh_xy[k][j].append(eval(self.lineEdit11[self.nx*self.ny*k+m].text()))
                    m+=1

        if num == (len(self.types)-1):
            self.close() 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window7(QWidget):
    def __init__(self,nxya,na,parent=None):
        super(Window7, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        self.setMinimumSize(QtCore.QSize(100, 100))
        self.nxya = nxya
        self.lab6  = [0]*nxya*nxya*na
        self.lineEdit10 = [0]*nxya*nxya*na
        self.Vect6 = [0]*nxya*nxya*na
        self.assembly = []
        self.layout = QVBoxLayout(self)
        typetab = QTabWidget(self)  
        typetab.setFont(QtGui.QFont("Sanserif", 10)) 
        self.types =  []
        for i in range(na):
            self.types.append("Assembly %s" %(i+1))
        num=0
        for name in self.types:
            tab =  QWidget()
            typetab.addTab(tab, name[0:10])
            typetablayout = QGridLayout(tab)
            m = 0
            for j in range(nxya):
                self.lab6[j] = QLabel("%s" %(j+1))
                self.lab6[j].setAlignment(Qt.AlignCenter)
                typetablayout.addWidget(self.lab6[j], j+7, 0)
                for i in range(nxya):
                    self.lab6[i] = QLabel("%s" %(i+1))
                    self.lab6[i].setAlignment(Qt.AlignCenter)
                    typetablayout.addWidget(self.lab6[i], 6, i+1)
                    self.lineEdit10[nxya*nxya*num+m] = QLineEdit()
                    typetablayout.addWidget(self.lineEdit10[nxya*nxya*num+m], j+7, i+1)
                    self.lineEdit10[nxya*nxya*num+m].insert(str(self.Vect6[nxya*nxya*num+m]))
                    m+=1
            self.layout.addWidget(typetab)
            self.setLayout(self.layout)
            #réer un Bouton
            if num == (len(self.types)-1):
                self.bouton = QPushButton(u"Save and Colse")
                self.bouton.clicked.connect(partial(self.save7,num))
                self.layout.addWidget(self.bouton)
            num+=1
    def save7(self,num):
	#connexion avec la fenetre main
        del self.assembly[:]
        for k in range(len(self.types)):
            m=0
            self.assembly.append([])
            for j in range(self.nxya):
                self.assembly[k].append([])
                for i in range(self.nxya):
                    self.assembly[k][j].append(eval(self.lineEdit10[self.nxya*self.nxya*k+m].text()))
                    m+=1
        if num == (len(self.types)-1):
            self.close() 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window8(QWidget):
    def __init__(self,na,parent=None):
        super(Window8, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        self.setMinimumSize(QtCore.QSize(100, 100))
        self.na = na
        self.lab4  = [0]*na*na
        self.lineEdit8 = [0]*na*na
        self.core   = [0]*na*na
        self.layout = QVBoxLayout(self)    
        typetablayout = QGridLayout()
        #Main group box
        box = QGroupBox()
        #self.main_group_box.setStyleSheet("QGroupBox{font-size: 10px}")
        box.setTitle("Core geometry:")
        box.setStyleSheet('QGroupBox:title {color: blue;}')
        #self.main_group_box.setLayout(self.layout)
        m = 0
        for j in range(na):
            self.lab4[j] = QLabel("%s" %(j+1))
            self.lab4[j].setAlignment(Qt.AlignCenter)
            typetablayout.addWidget(self.lab4[j], j+7, 0)
            for i in range(na):
                self.lab4[i] = QLabel("%s" %(i+1))
                self.lab4[i].setAlignment(Qt.AlignCenter)
                typetablayout.addWidget(self.lab4[i], 6, i+1)
                self.lineEdit8[m] = QLineEdit()
                typetablayout.addWidget(self.lineEdit8[m], j+7, i+1)
                self.lineEdit8[m].insert(str(self.core[m]))
                m+=1

        box.setLayout(typetablayout)
        self.layout.addWidget(box)
        self.setLayout(self.layout)

        #réer un Bouton
        self.bouton = QPushButton(u"Save and Colse")
        self.bouton.clicked.connect(self.save8)
        self.layout.addWidget(self.bouton)

    def save8(self):
	#connexion avec la fenetre main
        del self.core[:] 
        m=0
        for j in range(self.na):
            self.core.append([])
            for i in range(self.na):
                self.core[j].append(eval(self.lineEdit8[m].text()))   
                m+=1
        self.close()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# An Example application QWidget containing the textEdit_4 to redirect stdout to
class Application(QtWidgets.QMainWindow):
    from app.func import visualization1D,visualization2D,P_Ordinate,powerpf,power,plot
    def __init__(self, title= "Default", parent=None):
        super(Application, self).__init__(parent)
        sys.stdout = EmittingStream(textWritten=self.normalOutputWritten)
        sys.stderr = EmittingStream(textWritten=self.normalOutputWritten)
        self.title = title
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowTitle(self.title)
        self._initButtons()
        self.isFileOpen = False
        self.isFileCreate = False
        widget = QWidget()
        topFiller = QWidget()
        topFiller.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.infoLabel = QLabel(
                "<i>Choose a menu option, or right-click to invoke a context menu</i>",
                alignment=Qt.AlignCenter)
        self.infoLabel.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        bottomFiller = QWidget()
        bottomFiller.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        vbox = QVBoxLayout()
        vbox.setContentsMargins(5, 5, 5, 5)
        vbox.addWidget(topFiller)
        vbox.addWidget(self.infoLabel)
        vbox.addWidget(bottomFiller)
        widget.setLayout(vbox)
             
    def _initButtons(self):
        self.ui.action_Open.triggered.connect(self.openFile)
        self.ui.action_About.triggered.connect(self.about)
        self.ui.action_About_Qt.triggered.connect(self.aboutQt)
        self.ui.action_Help.triggered.connect(self.help)
        self.ui.action_Save.triggered.connect(self.saveFile)
        self.ui.action_Save_as.triggered.connect(self.save_asFile)
        self.ui.action_website.triggered.connect(self.webSite)
        self.ui.action_New.triggered.connect(self.NewFile)
        self.ui.action_Exit.triggered.connect(self.Exit)
        self.ui.action_New.triggered.connect(self.CloseFile)
        self.ui.action_Copy.triggered.connect(self.copy)
        self.ui.action_Paste.triggered.connect(self.paste)
        self.ui.action_Cut.triggered.connect(self.cut)
        self.ui.action_Run.triggered.connect(self.start_thread1)
        self.ui.pushButton_5.clicked.connect(self.start_thread1)
        self.ui.action_Compile_2.triggered.connect(self.compile)
        self.ui.action_Plot.triggered.connect(self.plot)
        self.ui.pushButton_22.clicked.connect(self.P_Ordinate)
        self.ui.pushButton_21.clicked.connect(self.plot)
        self.ui.pushButton_2.clicked.connect(self.visualization1D)
        self.ui.pushButton_12.clicked.connect(self.visualization2D)
        self.ui.pushButton_3.clicked.connect(self.plot)
        self.ui.pushButton_4.clicked.connect(self.compile)
        self.ui.pushButton_18.clicked.connect(self.compile)
        self.ui.pushButton_5.clicked.connect(self.run)
        self.ui.pushButton_19.clicked.connect(self.powerpf)
        self.ui.pushButton_6.clicked.connect(self.power)
        self.ui.pushButton_20.clicked.connect(self.run)
        self.ui.pushButton.clicked.connect(self.new1)
        self.ui.pushButton_7.clicked.connect(self.new2)
        self.ui.pushButton_26.clicked.connect(self.new5)
        self.ui.pushButton_23.clicked.connect(self.new6) 
        self.ui.pushButton_11.clicked.connect(self.data_up1D)
        self.ui.pushButton_16.clicked.connect(self.data_up2D)
        self.statusbar = self.statusBar()
        self.ui.textEdit_4.cursorPositionChanged.connect(self.cursorPosition_2)
        self.ui.textEdit_3.cursorPositionChanged.connect(self.cursorPosition)
        self.ui.radioButton_4.toggled.connect(lambda:self.methods_1D(self.ui.radioButton_4))
        self.ui.radioButton_2.toggled.connect(lambda:self.methods_1D(self.ui.radioButton_2))
        self.ui.radioButton_3.toggled.connect(lambda:self.methods_1D(self.ui.radioButton_3))  
        self.ui.radioButton.toggled.connect(lambda:self.methods_2D(self.ui.radioButton)) 
        self.ui.spinBox_7.setRange(-12, -1)
        self.ui.spinBox_7.setValue(-6)
        self.ui.spinBox_4.setSingleStep(2)
        self.ui.spinBox_12.setSingleStep(2)
        self.ui.spinBox_16.setRange(-12, -1)
        self.ui.spinBox_16.setValue(-6)
        self.ui.spinBox_6.setRange(1, 1000) 
        self.ui.spinBox_6.setValue(200)
        self.ui.spinBox_15.setRange(1, 1000) 
        self.ui.spinBox_15.setValue(200)
        self.ui.spinBox_18.setRange(1, 5000)

        # Select active method
        M00 = open('app/link/script00.py', "r" ).read() 
        if M00 == 'CP1D':
            self.ui.radioButton_2.setChecked(True) 
        elif M00 == 'SN1D':
            self.ui.radioButton_3.setChecked(True)
        elif M00 == 'MOC1D':
            self.ui.radioButton_4.setChecked(True)       
        elif M00 == 'SN2D':
            self.ui.radioButton.setChecked(True) 

        # Multi-group Cross Sections 1D
        self.ui.comboBox_13.currentIndexChanged.connect(self.Cross1D)
        M12 = open('app/link/script13.py', "r" ).read() 
        if M12 == 'TotalXS':
            self.ui.comboBox_13.setCurrentText('TotalXS') 
        elif M12 == 'FissionXS':
            self.ui.comboBox_13.setCurrentText('FissionXS')
        elif M12 == 'NuFissionXS':
            self.ui.comboBox_13.setCurrentText('NuFissionXS')
        elif M12 == 'ScatterMatrixXS':
            self.ui.comboBox_13.setCurrentText('ScatterMatrixXS')
        elif M12 == 'Chi':
            self.ui.comboBox_13.setCurrentText('ChiXS')

        # Multi-group Cross Sections 2D
        self.ui.comboBox_12.currentIndexChanged.connect(self.Cross2D)
        M12 = open('app/link/script11.py', "r" ).read() 
        if M12 == 'TotalXS':
            self.ui.comboBox_12.setCurrentText('TotalXS') 
        elif M12 == 'FissionXS':
            self.ui.comboBox_12.setCurrentText('FissionXS')
        elif M12 == 'NuFissionXS':
            self.ui.comboBox_12.setCurrentText('NuFissionXS')
        elif M12 == 'ScatterMatrixXS':
            self.ui.comboBox_12.setCurrentText('ScatterMatrixXS')
        elif M12 == 'Chi':
            self.ui.comboBox_12.setCurrentText('ChiXS')

        # Constructive geometry 1D
        self.ui.comboBox_6.currentIndexChanged.connect(self.ConstGeom1D)
        M10 = open('app/link/script12.py', "r" ).read()
        if M10 == 'Pin Cell':
            self.ui.comboBox_6.setCurrentText('Pin Cell') 
        elif M10 == 'Assembly':
            self.ui.comboBox_6.setCurrentText('Assembly')
        elif M10 == 'Core':
            self.ui.comboBox_6.setCurrentText('Core')

        # Constructive geometry 2D
        self.ui.comboBox_7.currentIndexChanged.connect(self.ConstGeom2D)
        M10 = open('app/link/script10.py', "r" ).read()
        if M10 == 'Pin Cell':
            self.ui.comboBox_7.setCurrentText('Pin Cell') 
        elif M10 == 'Assembly':
            self.ui.comboBox_7.setCurrentText('Assembly')
        elif M10 == 'Core':
            self.ui.comboBox_7.setCurrentText('Core')
        # Geometry
        self.ui.comboBox.currentIndexChanged.connect(self.geometry)
        M01 = open('app/link/script01.py', "r" ).read() 
        if M01 == 'Slab Geometry':
            self.ui.comboBox.setCurrentText('Slab Geometry') 
        elif M01 == 'Cylindrical Geometry':
            self.ui.comboBox.setCurrentText('Cylindrical Geometry')
        elif M01 == 'Spherical Geometry':
            self.ui.comboBox.setCurrentText('Spherical Geometry')
        # Boundary Conditions 1D
        self.ui.comboBox_2.currentIndexChanged.connect(self.bc)
        M02 = open('app/link/script02.py', "r" ).read()
        if M02 == 'Vacuum':
            self.ui.comboBox_2.setCurrentText('Vacuum') 
        elif M02 == 'Vacuum Reflective':
            self.ui.comboBox_2.setCurrentText('Vacuum Reflective')
        elif M02 == 'Reflective Vacuum':
            self.ui.comboBox_2.setCurrentText('Reflective Vacuum')
        else:
            self.ui.comboBox_2.setCurrentText('Reflective')
        # Boundary Conditions 2D
        self.ui.comboBox_5.currentIndexChanged.connect(self.bc2d1)
        M00 = open('app/link/script06.py', "r" ).read()
        if M00 == 'Vacuum Right':
            self.ui.comboBox_5.setCurrentText('Vacuum Right') 
        else:
            self.ui.comboBox_5.setCurrentText('Reflective Right')

        self.ui.comboBox_9.currentIndexChanged.connect(self.bc2d2)
        M00 = open('app/link/script07.py', "r" ).read()
        if M00 == 'Reflective Left':
            self.ui.comboBox_9.setCurrentText('Reflective Left') 
        else:
            self.ui.comboBox_9.setCurrentText('Vacuum Left')

        self.ui.comboBox_11.currentIndexChanged.connect(self.bc2d3)
        M00 = open('app/link/script08.py', "r" ).read()
        if M00 == 'Vacuum Bottom':
            self.ui.comboBox_11.setCurrentText('Vacuum Bottom') 
        else:
            self.ui.comboBox_11.setCurrentText('Reflective Bottom')

        self.ui.comboBox_10.currentIndexChanged.connect(self.bc2d4)
        M00 = open('app/link/script09.py', "r" ).read()
        if M00 == 'Vacuum Top':
            self.ui.comboBox_10.setCurrentText('Vacuum Top') 
        else:
            self.ui.comboBox_10.setCurrentText('Reflective Top')
        # Approximation Scheme 1
        self.ui.comboBox_3.currentIndexChanged.connect(self.asc1)
        self.ui.comboBox_3.setToolTip('Choose Discritization Scheme for the Discrete Ordinates <b>(SN)</b> method')
        M03 = open('app/link/script03.py', "r" ).read()
        if M03 == 'Step Difference':
            self.ui.comboBox_3.setCurrentText('Step Difference') 
        elif M03 == 'Diamond Difference':
            self.ui.comboBox_3.setCurrentText('Diamond Difference')
        # Approximation Scheme 2
        self.ui.comboBox_4.currentIndexChanged.connect(self.asc2)
        self.ui.comboBox_4.setToolTip('Choose Discritization Scheme for the Method Of Characteristics <b>(MOC)</b>')
        M04 = open('app/link/script04.py', "r" ).read()
        if M04 == 'Step Characteristics':
            self.ui.comboBox_4.setCurrentText('Step Characteristics') 
        elif M04 == 'DD0':
            self.ui.comboBox_4.setCurrentText('DD0')
        else:
            self.ui.comboBox_4.setCurrentText('DD1')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def new1(self):
        npc     = self.ui.spinBox_3.value() 
        ngauss  = self.ui.spinBox_4.value()
        order   = self.ui.spinBox_5.value()
        nxpc    = self.ui.spinBox_22.value()
        npca    = self.ui.spinBox_21.value()
        na      = self.ui.spinBox_20.value()
        apc     = self.ui.spinBox_19.value()

        M05 = open('app/link/script12.py', "r" ).read() 
        if M05 == 'Pin Cell':
            if  npc  == 0:
                QMessageBox.warning(self, "Warning", "Enter the pin cells number")
            elif  nxpc == 0:
                QMessageBox.warning(self, "Warning", "Enter the x mesh pin cell number \n(Each pin cell is approximated by a X cartesian grid)")
            elif npc != 0 and nxpc != 0:
                self.wind1 = Window1(npc,nxpc)
                self.wind1.show()
        elif M05 == 'Assembly':
            if  na  == 0:
                QMessageBox.warning(self, "Warning", "Enter the assemblies number")
            elif  npca  == 0:
                QMessageBox.warning(self, "Warning", "Enter the number of pin cell in each assembly \n(Each assembly contains a set of pin cells)")
            elif na != 0 and npca !=0:
                self.wind2 = Window2(npca,na)
                self.wind2.show()
                
        elif M05 == 'Core':
            if na == 0:
                QMessageBox.warning(self, "Warning", "Enter the assemblies number")
            elif na != 0:
                self.wind3 = Window3(na)
                self.wind3.show()
        else:
            pass

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function for Multi-group Cross Sections 1D
    def new2(self):
        M06 = open('app/link/script13.py', "r" ).read() 
        ngroup  = self.ui.spinBox.value()
        nmat    = self.ui.spinBox_2.value()
        order   = self.ui.spinBox_5.value()
        if  ngroup == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of energy group")
        elif nmat == 0:
            QMessageBox.warning(self, "Warning", "Enter the materials number")
        else:          
            if M06 == 'TotalXS':
                self.wind9 = Window9(ngroup,nmat)
                self.wind9.show()
            elif M06 == 'FissionXS':
                self.wind10 = Window10(ngroup,nmat)
                self.wind10.show()
            elif M06 == 'NuFissionXS':
                self.wind11 = Window11(ngroup,nmat)
                self.wind11.show()
            elif M06 == 'ScatterMatrixXS':
                self.wind12 = Window12(ngroup,nmat) 
                self.wind12.show()
            elif M06 == 'ChiXS':
                self.wind13 = Window13(ngroup,nmat)
                self.wind13.show()
            else:
                print ('error')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def new5(self):
        np   = self.ui.spinBox_9.value()
        na   = self.ui.spinBox_11.value()
        nmat = self.ui.spinBox_10.value()
        nxy  = self.ui.spinBox_14.value()
        nxya   = self.ui.spinBox_17.value()
        M05 = open('app/link/script10.py', "r" ).read() 
        if M05 == 'Pin Cell':
            if nmat == 0:
                QMessageBox.warning(self, "Warning", "Enter the materials number")
            elif  np  == 0:
                QMessageBox.warning(self, "Warning", "Enter the pin cells number")
            elif  nxy == 0:
                QMessageBox.warning(self, "Warning", "Enter the x-y mesh pin cell number \n(Each pin cell is approximated by a X x Y cartesian grid)")
            elif np != 0 and nxy != 0 and nmat != 0:
                self.wind6 = Window6(nxy,nxy,np)
                self.wind6.show()
        elif M05 == 'Assembly':
            if  na  == 0:
                QMessageBox.warning(self, "Warning", "Enter the assemblies number")
            elif  nxya  == 0:
                QMessageBox.warning(self, "Warning", "Enter the x-y mesh assembly number \n(Each assembly contains a set of pin cells)")
            elif na != 0 and nxya !=0:
                self.wind7 = Window7(nxya,na)
                self.wind7.show()
                
        elif M05 == 'Core':
            if na == 0:
                QMessageBox.warning(self, "Warning", "Enter the assemblies number")
            elif na != 0:
                self.wind8 = Window8(na)
                self.wind8.show()
        else:
            pass
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function for Multi-group Cross Sections 2D
    def new6(self):
        M06 = open('app/link/script11.py', "r" ).read() 
        ngroup   = self.ui.spinBox_8.value()
        nmat     = self.ui.spinBox_10.value()
        if  ngroup == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of energy group")
        elif nmat == 0:
            QMessageBox.warning(self, "Warning", "Enter the materials number")
        else:          
            if M06 == 'TotalXS':
                self.wind9 = Window9(ngroup,nmat)
                self.wind9.show()
            elif M06 == 'FissionXS':
                self.wind10 = Window10(ngroup,nmat)
                self.wind10.show()
            elif M06 == 'NuFissionXS':
                self.wind11 = Window11(ngroup,nmat)
                self.wind11.show()
            elif M06 == 'ScatterMatrixXS':
                self.wind12 = Window12(ngroup,nmat)
                self.wind12.show()
            elif M06 == 'ChiXS':
                self.wind13 = Window13(ngroup,nmat)
                self.wind13.show()
            else:
                print ('error')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def data_up1D(self):
        Geometry_type = open('app/link/script01.py', "r" ).read()
        ngroup  = self.ui.spinBox.value()
        nmat    = self.ui.spinBox_2.value()
        order   = self.ui.spinBox_5.value()
        npc     = self.ui.spinBox_3.value() 
        ngauss  = self.ui.spinBox_4.value()
        order   = self.ui.spinBox_5.value()
        nxpc    = self.ui.spinBox_22.value()
        npca    = self.ui.spinBox_21.value()
        na      = self.ui.spinBox_20.value()
        apc     = self.ui.spinBox_19.value()
        max_it  = self.ui.spinBox_6.value()
        Tol     = self.ui.spinBox_7.value()

        bc = self.ui.comboBox_2.currentText() 

        if  ngroup == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres Energy Group Number")
        elif nmat == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres Number of Materials")
        elif npc == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres pin cells number")
        elif nxpc == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres x mesh pin cell number") 
        elif npca == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres x mesh assembly number") 
        elif na == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres assemblies number")
        elif ngauss == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres SN Order")
        elif apc== 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres active pin cells number")
        else: 
            wind1 = Window1(npc,nxpc)
            wind2 = Window2(npca,na)
            wind3 = Window3(na)
            wind9 = Window9(ngroup,nmat)
            wind10 = Window10(ngroup,nmat)
            wind11 = Window11(ngroup,nmat)
            wind12 = Window12(ngroup,nmat)
            wind13 = Window13(ngroup,nmat)

            filename = open("app/input/input.json",'w')
            open('app/link/script.dir', "w" ).write(os.path.abspath(os.path.dirname( __file__)) +'/input/input.json')
            filename.write('{ \n  "data": { \n    "parameter": { \n      "id": 100,')
            filename.write('\n      "Total number of energy groups": '+ str(ngroup) +',')
            filename.write('\n      "Total number of Materials": '+ str(nmat) +',')
            filename.write('\n      "Total number of pin cells": '+ str(npc) +',')
            filename.write('\n      "Total number of assemblies": '+ str(na) +',')
            filename.write('\n      "Core": '+ str(self.wind3.core) +',')
            filename.write('\n      "Total number of active pin cells": '+ str(apc) +',')
            filename.write('\n      "Number of angular discretizations": '+ str(ngauss) +',')
            filename.write('\n      "The l-order Legendre polynomial": '+ str(order) +',')
            filename.write('\n      "Boundary conditions": '+'"'+ str(bc)+'"' +',')
            filename.write('\n      "Maximum number of iterations": '+ str(max_it) +',')
            filename.write('\n      "Criterion of Keff convergence": '+ '1.0E'+str(Tol))
            filename.write('\n    }, \n    "Assemblies": [')
            # Boucle sur assemblage
            for i in range(na):
                filename.write('\n      { \n        "id": '+ str(i+1) +', \n        "name": "Assembly '+ str(i+1) +'",') 
                filename.write('\n        "assembly": ' + str(self.wind2.assembly[i][:])) 
                if i == na-1:
                    filename.write('\n      }')
                else:
                    filename.write('\n      },')  
            filename.write('\n], \n    "PinCells": [')
            # Boucle sur les pin cell
            for i in range(npc):
                filename.write('\n      { \n        "id": '+ str(i+1) +', \n        "name": "Pin Cell '+ str(i+1) +'",') 

                if Geometry_type in ['Cylindrical Geometry','Spherical Geometry']:
                    filename.write('\n        "ray": ' + str(self.wind1.dx[i][:]) +',')
                else:
                    filename.write('\n        "width_x": ' + str(self.wind1.dx[i][:]) +',')
                filename.write('\n        "mat_fill": '+ str(self.wind1.regmat[i][:]) +',')
                filename.write('\n        "fine_mesh": '+ str(self.wind1.nfmesh_x[i][:]))
                if i == npc-1:
                    filename.write('\n      }')
                else:
                    filename.write('\n      },')
            filename.write('\n], \n    "materials": [')
            # Boucle sur les matériaux
            for i in range(nmat):
                filename.write('\n      { \n        "id": '+ str(i+1) +', \n        "name": "material '+ str(i+1) +'",') 
                filename.write('\n        "XSTotal": ' + str(self.wind9.SigT[i][:]) +','
                                   + '\n        "XSFission": '+ str(self.wind10.SigF[i][:])+','
                                   + '\n        "XSNuFission": '+ str(self.wind11.NuSigF[i][:])+',')
                filename.write('\n        "XSScatter Matrix":'+str(self.wind12.SigS[i][:])+','+ 
                 '\n        "XSChi":  '+str(self.wind13.Chi[i][:]))
                if i == nmat-1:
                    filename.write('\n      }')
                else:
                    filename.write('\n      },')  
            # Fin Boucle
            filename.write('\n    ]  \n  }  \n}') 
            filename.close() 
            fh = open("app/input/input.json","r")  
            self.ui.textEdit_4.setText(fh.read())    
            fh.close()  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def data_up2D(self):   
        bc = [0]*4
        ngroup   = self.ui.spinBox_8.value()
        np       = self.ui.spinBox_9.value()   
        nmat     = self.ui.spinBox_10.value() 
        na       = self.ui.spinBox_11.value() 
        nxy      = self.ui.spinBox_14.value()
        nxya     = self.ui.spinBox_17.value()
        ordinate = self.ui.spinBox_12.value()
        apc      = self.ui.spinBox_18.value()
        order    = self.ui.spinBox_13.value() 
        max_it   = self.ui.spinBox_15.value()
        Tol      = self.ui.spinBox_16.value()
        if self.ui.comboBox_5.currentText()== 'Vacuum Right':
            bc[0]= "VR"
        else:
            bc[0]= "RR"
        if self.ui.comboBox_9.currentText()== 'Vacuum Left':
            bc[1]= "VL"
        else:
            bc[1]= "RL"
        if self.ui.comboBox_10.currentText()== 'Vacuum Bottom':
            bc[2]= "VB"
        else:
            bc[2]= "RB"
        if self.ui.comboBox_11.currentText()== 'Vacuum Top':
            bc[3]= "VT"
        else:
            bc[3]= "RT"
        if  ngroup == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres Energy Group Number")
        elif nmat == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres Number of Materials")
        elif np == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres pin cells number")
        elif nxy == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres x-y mesh pin cell number") 
        elif nxya == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres x-y mesh assembly number") 
        elif na == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres assemblies number")
        elif ordinate == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres SN Order")
        elif apc== 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres active pin cells number")
        else: 
            wind6  = Window6(nxy,nxy,np)
            wind7  = Window7(nxya,na)
            wind8  = Window8(na)
            wind9  = Window9(ngroup,nmat) 
            wind10 = Window10(ngroup,nmat)
            wind11 = Window11(ngroup,nmat)
            wind12 = Window12(ngroup,nmat)
            wind13 = Window13(ngroup,nmat)
            filename = open("app/input/input.json",'w')
            open('app/link/script.dir', "w" ).write(os.path.abspath(os.path.dirname( __file__)) +'/input/input.json')
            filename.write('{ \n  "data": { \n    "parameter": { \n      "id": 100,')
            filename.write('\n      "Total number of energy groups": '+ str(ngroup) +',')
            filename.write('\n      "Total number of materials": '+ str(nmat) +',')
            filename.write('\n      "Total number of pin cells": '+ str(np) +',')
            filename.write('\n      "Total number of assemblies": '+ str(na) +',')
            filename.write('\n      "Core": '+ str(self.wind8.core) +',')
            filename.write('\n      "Total number of active pin cells": '+ str(apc) +',')
            filename.write('\n      "Number of angular discretizations": '+ str(ordinate) +',')
            filename.write('\n      "The l-order Legendre polynomial": '+ str(order) +',')
            filename.write('\n      "Boundary conditions": '+ '['+'"'+bc[0]+'"'+','+'"'+bc[1]+'"'+','+'"'+bc[2]+'"'+','+'"'+bc[3]+'"'+']' +',')
            filename.write('\n      "Maximum number of iterations": '+ str(max_it) +',')
            filename.write('\n      "Criterion of Keff convergence": '+ '1.0E'+str(Tol))
            filename.write('\n    }, \n    "Assemblies": [')
            # Boucle sur assemblage
            for i in range(na):
                filename.write('\n      { \n        "id": '+ str(i+1) +', \n        "nom": "Assembly '+ str(i+1) +'",') 
                filename.write('\n        "assembly": ' + str(self.wind7.assembly[i][:])) 
                if i == na-1:
                    filename.write('\n      }')
                else:
                    filename.write('\n      },')  
            filename.write('\n], \n    "PinCells": [')
            # Boucle sur les pin cell
            for i in range(np):
                filename.write('\n      { \n        "id": '+ str(i+1) +', \n        "nom": "Pin Cell '+ str(i+1) +'",') 
                filename.write('\n        "width_x": ' + str(self.wind6.dx[i][:]) +','
                                   + '\n        "width_y": '+ str(self.wind6.dy[i][:])+',')
                filename.write('\n        "mat_fill": '+ str(self.wind6.regmat[i][:]) +',')
                filename.write('\n        "fine_mesh": '+ str(self.wind6.nfmesh_xy[i][:]))
                if i == np-1:
                    filename.write('\n      }')
                else:
                    filename.write('\n      },')
            filename.write('\n], \n    "materials": [')
            # Boucle sur les matériaux
            for i in range(nmat):
                filename.write('\n      { \n        "id": '+ str(i+1) +', \n        "nom": "material '+ str(i+1) +'",') 
                filename.write('\n        "XSTotal": ' + str(self.wind9.SigT[i][:]) +','
                                   + '\n        "XSFission": '+ str(self.wind10.SigF[i][:])+','
                                   + '\n        "XSNuFission": '+ str(self.wind11.NuSigF[i][:])+',')
                filename.write('\n        "XSScatter Matrix":'+str(self.wind12.SigS[i][:])+','+ 
                 '\n        "XSChi":  '+str(self.wind13.Chi[i][:]))
                if i == nmat-1:
                    filename.write('\n      }')
                else:
                    filename.write('\n      },')  
            # Fin Boucle
            filename.write('\n    ]  \n  }  \n}') 
            filename.close() 
            fh = open("app/input/input.json","r")  
            self.ui.textEdit_4.setText(fh.read())    
            fh.close()  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def aboutQt(self):
        self.infoLabel.setText("Invoked <b>Help|About Qt</b>") 

    def Cross1D(self,b):
        open('app/link/script13.py', "w" ).write(self.ui.comboBox_13.currentText())
 
    def Cross2D(self,b):
        open('app/link/script11.py', "w" ).write(self.ui.comboBox_12.currentText()) 

    def ConstGeom2D(self,b):
        open('app/link/script10.py', "w" ).write(self.ui.comboBox_7.currentText()) 

    def ConstGeom1D(self,b):
        open('app/link/script12.py', "w" ).write(self.ui.comboBox_6.currentText())

    def geometry(self,b):
        open('app/link/script01.py', "w" ).write(self.ui.comboBox.currentText()) 
	
    def bc(self,b):	
        open('app/link/script02.py', "w" ).write(self.ui.comboBox_2.currentText()) 

    def bc2d1(self,b):	
        open('app/link/script06.py', "w" ).write(self.ui.comboBox_5.currentText())

    def bc2d2(self,b):	
        open('app/link/script07.py', "w" ).write(self.ui.comboBox_9.currentText())

    def bc2d3(self,b):	
        open('app/link/script08.py', "w" ).write(self.ui.comboBox_11.currentText())

    def bc2d4(self,b):	
        open('app/link/script09.py', "w" ).write(self.ui.comboBox_10.currentText())

    def asc1(self,b):	
        open('app/link/script03.py', "w" ).write(self.ui.comboBox_3.currentText()) 

    def asc2(self,b):	
        open('app/link/script04.py', "w" ).write(self.ui.comboBox_4.currentText()) 

    def methods_1D(self,b):
        if b.isChecked() == True: 
            self.ui.radioButton.setChecked(False)
            open('app/link/script00.py', "w" ).write(b.text())
        else:
            open('app/link/script00.py', "w" ).close()
    def methods_2D(self,b):
        if b.isChecked() == True:
            self.ui.radioButton_2.setChecked(False)
            self.ui.radioButton_3.setChecked(False)
            self.ui.radioButton_4.setChecked(False)
            open('app/link/script00.py', "w" ).write(b.text()) 
        else:
            open('app/link/script00.py', "w" ).close()

    def cursorPosition(self):
        cursor = self.ui.textEdit_3.textCursor()

        # Mortals like 1-indexed things
        line = cursor.blockNumber() + 1
        col = cursor.columnNumber()
        self.statusbar.showMessage("Line: {} | Column: {} --> Output Window".format(line,col))

    def cursorPosition_2(self):
        cursor = self.ui.textEdit_4.textCursor()
        # Mortals like 1-indexed things
        line = cursor.blockNumber() + 1
        col = cursor.columnNumber()
        self.statusbar.showMessage("Line: {} | Column: {} --> Input Window".format(line,col))

    def copy(self):
        cursor=self.ui.textEdit_4.textCursor()
        textSelected = cursor.selectedText()
        self.copiedtext=textSelected
    
    def paste(self):
        self.ui.textEdit_4.append(self.copiedtext)

    def cut(self):
        cursor=self.ui.textEdit_4.textCursor()
        textSelected=cursor.selectedText()
        self.copiedtext=textSelected
        self.ui.textEdit_4.cut()

    def NewFile(self):
        self.isFileCreate = True
        self.CloseFile()
        self.ui.textEdit_4.show()

    def CloseFile(self):
        if self.ui.textEdit_4.toPlainText() == "":
            pass
        else:
            reply = QMessageBox.question(self,"", "Are you sure you want to close this file ?", QMessageBox.Yes | QMessageBox.No)
            if reply == QMessageBox.Yes:
                self.ui.textEdit_4.clear()
                self.setWindowTitle(self.title)
                self.isFileCreate = True
                self.isFileOpen = False
            else:
                pass

    def Exit(self):
        """Generate 'question' dialog on clicking 'X' button in title bar.
        Reimplement the closeEvent() event handler to include a 'Question'
        dialog with options on how to proceed - Save, Close, Cancel buttons
        """
        reply = QMessageBox.question(self, "Message",
            "Are you sure you want to quit ?", QMessageBox.Yes, QMessageBox.No)

        if reply == QMessageBox.Yes:
            qapp.quit()
        else:
            pass

    def openFile(self):
        self.filename = QtWidgets.QFileDialog.getOpenFileName(self,'Open File', "~", "*.json")[0]
        open('app/link/script.dir', "w" ).write(str(self.filename))
        if self.filename != "":
            self.isFileOpen = True
            file = open(self.filename,"r")
            self.setWindowTitle(self.title + ':' + self.filename)
            self.ui.textEdit_4.show()
            self.ui.textEdit_4.setText(file.read())
            file.close()

    def saveFile(self):
        if self.isFileOpen == False or self.isFileCreate == False:
            pass

        if self.isFileOpen == True:
            file = open(self.filename,"w")
            self.setWindowTitle(self.title + ':' + self.filename)
            file.write(self.ui.textEdit_4.toPlainText())
            file.close()
            self.ui.statusbar.showMessage(self.parseFileName() + " has been saved " + " at " +
                                          time.strftime('%d/%m/%y %H:%M', time.localtime()), 4600)
        elif self.isFileCreate == True:
            self.isFileOpen = True
            self.filename = QtWidgets.QFileDialog.getSaveFileName(self,'Save File', "data.json", "*.json")[0]
            if self.filename != "":
                file = open(self.filename,"w")
                self.setWindowTitle(self.title + ':' + self.filename)
                file.write(self.ui.textEdit_4.toPlainText())
                file.close()
                self.ui.statusbar.showMessage(self.parseFileName()  + " has been saved " + " at " +
                                              time.strftime('%d/%m/%y %H:%M', time.localtime()), 4600)
        else:
            pass

    def save_asFile(self):
        self.filename = QtWidgets.QFileDialog.getSaveFileName(self,'Save as File', "data.json", "*.json")[0]
        if self.filename != "":
            file = open(self.filename,"w")
            self.setWindowTitle(self.title + ':' + self.filename)
            file.write(self.ui.textEdit_4.toPlainText())
            file.close()
            self.ui.statusbar.showMessage(self.parseFileName()  + " has been saved " + " at " +
                                          time.strftime('%d/%m/%y %H:%M', time.localtime()), 4600)
        else:
            pass

    def parseFileName(self):
        filename = self.filename.split("/")
        self.fname = filename[-1]
        return self.fname        

    def about(self):
        QMessageBox.about(self, "About OpenNTP",
                "<center> Python GUI Programming Using Qt <center>\n \n" 
                "<center> This project was developed by <center>\n"
                "<center> Mohamed LAHDOUR & Tarek EL Bardouni. <center>\n"
                "<center> Departement of Physics, Laboratory of Radiation"
                                  " & Nuclear Systems <center>\n"
                "<center> University Abdelmalek Essaadi, Faculty of sciences"
                                  "Tetouan (Morocco). <center>")

    def help(self):
        QMessageBox.about(self, "OpenNTP",
                "<center> For Help Contact Us: <center>\n \n \n" 
                "<center> mlahdour@uae.ac.ma  &  tarekbardouni@uae.ma <center>")

    def webSite(self):
        url ="https://github.com/mohamedlahdour/NTP-ERSN/releases"
        self.ui.statusbar.showMessage('Loading url...', 5600)
        webbrowser.open(url)

    @pyqtSlot()
    def start_thread1(self):
        self.thread = QThread()
        self.long_running_thing = LongRunningThing()
        self.long_running_thing.moveToThread(self.thread)
        self.thread.started.connect(self.long_running_thing.run)
        self.thread.start()
        self.thread.quit()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    @pyqtSlot()
    def compile(self):
        """exécuté lors du clic sur le bouton
        """
        #print ('Connecting process')
        self.process = QtCore.QProcess(self)
        self.process.readyReadStandardOutput.connect(self.stdoutReady)
        self.process.readyReadStandardError.connect(self.stderrReady)
        #print ('Starting process')
        self.process.start('python3', ['app/compile.py'])

    @pyqtSlot()
    def run(self):
        #print ('Connecting process')
        self.process = QtCore.QProcess(self)
        self.process.readyReadStandardOutput.connect(self.stdoutReady)
        self.process.readyReadStandardError.connect(self.stderrReady)
        #print ('Starting process')
        self.process.start('python3', ['main.py'])
        
    def append(self, text):
        cursor = self.ui.textEdit_3.textCursor()
        cursor.insertText(text)
        self.ui.textEdit_3.ensureCursorVisible()

    def append_text(self,text):
        cursor = self.ui.textEdit_3.textCursor()
        cursor.insertText(text)
        self.ui.textEdit_3.ensureCursorVisible()

    def stdoutReady(self):
        text = bytearray(self.process.readAllStandardOutput())
        text = text.decode("ascii")
        self.append_text(text)


    def stderrReady(self):
        text = bytearray(self.process.readAllStandardError())
        text = text.decode("ascii")
        self.append_text(text)

    def normalOutputWritten(self,text):
        cursor = self.ui.textEdit_3.textCursor()
        cursor.insertText(text)
        
if __name__ == '__main__': 
    qapp = QApplication(sys.argv) 
    app = Application(u'OpenNTP')
    app.show()
    qapp.exec_()
        



