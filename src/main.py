# -*- coding: utf-8 -*-
import sys
import os
from multiprocessing import Queue
from PyQt5 import QtCore, QtGui, QtWidgets
from mainwindow import Ui_MainWindow
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
"""
=======
Buttons
=======

Constructing a simple button GUI to modify a sine wave.

The ``next`` and ``previous`` button widget helps visualize the wave with
new frequencies.
"""

class Index(object):
    def __init__(self,ng,ax,clb,parent=None):
        self.ng = ng
        self.ind = ng
        self.freqs = np.arange(0, ng, 1)
        self.ax  = ax
        self.clb = clb
    def next(self, event):
        self.ind += 1
        i = self.ind % len(self.freqs)
        if self.ind > self.ng:
            pass
        else:
            if self.ind <= 1:
                self.ind = 2
                if self.ng ==1:
                    self.ind = 1
            else:
                pass
            data = np.loadtxt('app/Output/FG'+str(abs(self.ind)))
            max_columns = len(data[0]) - 1
            max_rows = len(data)
            x = [data[rownum+1][0] for rownum in range(max_rows-1)] 
            y = [data[0][colnum + 1] for colnum in range(max_columns)]
            z = [[data[rownum+1][colnum + 1] for rownum in range(max_rows-1)] for colnum in range(max_columns)]
            x = np.array(x)
            y = np.array(y)
            z = np.array(z)
            x,y = np.meshgrid(x, y, sparse=True)
            self.ax.set_title('Energy Group '+str(abs(self.ind)))
            im = self.ax.imshow(z, interpolation='bilinear', cmap='jet', 
            origin='lower', extent=[0, abs(x).max(), 0,abs(y).max()])
            self.clb.on_mappable_changed(im)


    def prev(self, event):
        self.ind -= 1
        i = self.ind % len(self.freqs)
        if self.ind < 1:
            pass
        else:
            if self.ind >= self.ng:
                self.ind = self.ng-1
                if self.ind == 0:
                    self.ind = 1
            else:
                pass
            data = np.loadtxt('app/Output/FG'+str(abs(self.ind)))
            max_columns = len(data[0]) - 1
            max_rows = len(data)
            x = [data[rownum+1][0] for rownum in range(max_rows-1)] 
            y = [data[0][colnum + 1] for colnum in range(max_columns)]
            z = [[data[rownum+1][colnum + 1] for rownum in range(max_rows-1)] for colnum in range(max_columns)]
            x = np.array(x)
            y = np.array(y)
            z = np.array(z)
            x,y = np.meshgrid(x, y, sparse=True)
            self.ax.set_title('Energy Group '+str(abs(self.ind)))
            im = self.ax.imshow(z, interpolation='bilinear', cmap='jet', 
            origin='lower', extent=[0, abs(x).max(), 0,abs(y).max()])
            self.clb.on_mappable_changed(im) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window1(QWidget):
    def __init__(self,nregion,nmat,parent=None):
        super(Window1, self).__init__(parent)
        self.setWindowTitle(u"Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        self.setMinimumSize(QtCore.QSize(100, 100))
        self.nr = nregion
        self.nm = nmat
        self.Delta  = [0]*nregion
        self.REGMAT = [0]*nregion
        self.NFMR   = [0]*nregion 
        self.lineEdit1 = [0]*100;self.lineEdit2 = [0]*100;self.lineEdit3 = [0]*100
        lab1  = [0]*100; lab2  = [0]*100; lab3  = [0]*100
        # créer un bouton
        self.bouton = QPushButton(u"Save")
        self.bouton.clicked.connect(self.save1)
        layout = QGridLayout()
        #------------------------------------
        Geometry_type = open('app/link/script01.py', "r" ).read()
        if Geometry_type in ['Cylindrical Geometry','Spherical Geometry']:
            tex= "<font color=blue > Ray for each Region per [cm]:</font>"
            tex1='Ray'
        else:
            tex= "<font color=blue > Size for each Region per [cm]:</font>"
            tex1='Size'
        topLabel1 = QLabel(tex) 
        lab1[0] = QLabel(tex1)
        for i in range(nregion):
            lab1[i+1] = QLabel("Region %s" %(i+1))
            self.lineEdit1[i] = QLineEdit()
            layout.addWidget(lab1[i+1], 1, i+1)
            layout.addWidget(self.lineEdit1[i], 2, i+1)
            self.lineEdit1[i].insert(str(self.Delta[i]))
        layout.addWidget(topLabel1,0,0)
        layout.addWidget(lab1[0], 2, 0)
        #------------------------------------
        topLabel2 = QLabel("<font color=blue > Number of fine meshes per region:</font>")    
        lab2[0] = QLabel('NFMR')
        for i in range(nregion):
            lab2[i+1] = QLabel("Region %s" %(i+1))
            self.lineEdit2[i] = QLineEdit()
            layout.addWidget(lab2[i+1], 4, i+1)
            layout.addWidget(self.lineEdit2[i], 5, i+1)
            self.lineEdit2[i].insert(str(self.REGMAT[i]))
        layout.addWidget(topLabel2,3,0)
        layout.addWidget(lab2[0], 5, 0)
        #------------------------------------
        topLabel3 = QLabel("<font color=blue > Which material fills each region:</font>")    
        lab3[0] = QLabel('Materials')
        for i in range(nregion):
            lab3[i+1] = QLabel("Region %s" %(i+1))
            self.lineEdit3[i] = QLineEdit()
            layout.addWidget(lab3[i+1], 7, i+1)
            layout.addWidget(self.lineEdit3[i], 8, i+1)
            self.lineEdit3[i].insert(str(self.NFMR[i]))
        layout.addWidget(topLabel3,6,0)
        layout.addWidget(lab3[0], 8, 0)
        layout.addWidget(self.bouton, 9, 0)
        self.setLayout(layout)

    def save1(self):
	#connexion avec la fenetre main
        del self.Delta[:]
        del self.REGMAT[:]
        del self.NFMR[:] 
        for i in range(self.nr):
            self.Delta.append(eval(self.lineEdit1[i].text()))
            self.NFMR.append(eval(self.lineEdit2[i].text()))
            self.REGMAT.append(eval(self.lineEdit3[i].text()))
        self.close()  

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window2(QWidget):
    def __init__(self,ngroup,nmat,parent=None):
        super(Window2, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        self.setMinimumSize(QtCore.QSize(100, 100))
        self.ng = ngroup
        self.nm = nmat
        self.Vect2  = [0]*ngroup*nmat
        self.SigT  = [0]*ngroup*nmat
        self.tab2  = [[0]*ngroup]*nmat
        self.lineEdit4 = [0]*ngroup*nmat 
        lab4  = [0]*ngroup*nmat
        layout = QGridLayout()
        #------------------------------------
        topLabel1 = QLabel("<font color=blue > Total Cross Section (SigT):</font>") 
        #if np.size(self.tab2) > np.size(self.SigT): 
        #    self.Vect2  = [0]*ngroup*nmat
        m = 0
        layout.addWidget(topLabel1,0,0)
        for j in range(nmat):
            lab4[j] = QLabel("Material %s" %(j+1))
            layout.addWidget(lab4[j], j+2, 0)
            for i in range(ngroup):
                lab4[i] = QLabel("G %s" %(i+1))
                layout.addWidget(lab4[i], 1, i+1)
                self.lineEdit4[m] = QLineEdit()
                layout.addWidget(self.lineEdit4[m], j+2, i+1)
                self.lineEdit4[m].insert(str(self.Vect2[m]))
                m=m+1
        #créer un bouton
        self.bouton = QPushButton(u"Save")
        self.bouton.clicked.connect(self.save2)
        layout.addWidget(self.bouton, j+3, 0)
        self.setLayout(layout)

    def save2(self):
        m=0
        self.SigT = []
        for i in range(self.nm):
            self.SigT.append([])
            for j in range(self.ng):
                self.SigT[i].append(eval(self.lineEdit4[m].text()))
                m+=1
        m = 0
        for i in range(self.nm):
            for j in range(self.ng):
                self.Vect2[m] = eval(self.lineEdit4[m].text())
                m+=1
        self.close()  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window3(QWidget):
    def __init__(self,ngroup,nmat,parent=None):
        super(Window3, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        self.setMinimumSize(QtCore.QSize(100, 100))
        self.ng = ngroup
        self.nm = nmat
        self.Vect3  = [0]*ngroup*nmat
        self.NuSigF  = [0]*ngroup*nmat
        self.tab3  = [[0]*ngroup]*nmat
        self.lineEdit5 = [0]*ngroup*nmat
        lab4  = [0]*ngroup*nmat
        layout = QGridLayout()
        #------------------------------------
        topLabel1 = QLabel("<font color=blue > NuFission Cross Section (NuSigF):</font>") 
        m = 0
        layout.addWidget(topLabel1,0,0)
        for j in range(nmat):
            lab4[j] = QLabel("Material %s" %(j+1))
            layout.addWidget(lab4[j], j+2, 0)
            for i in range(ngroup):
                lab4[i] = QLabel("G %s" %(i+1))
                layout.addWidget(lab4[i], 1, i+1)
                self.lineEdit5[m] = QLineEdit()
                layout.addWidget(self.lineEdit5[m], j+2, i+1)
                self.lineEdit5[m].insert(str(self.Vect3[m]))
                m=m+1
         # créer un bouton
        self.bouton = QPushButton(u"Save")
        self.bouton.clicked.connect(self.save3)
        layout.addWidget(self.bouton, j+3, 0)
        self.setLayout(layout)

    def save3(self):
        m=0
        self.Vect3 = [0]*self.nm*self.ng
        del self.NuSigF[:]
        for i in range(self.nm):
            self.NuSigF.append([])
            for j in range(self.ng):
                self.NuSigF[i].append(eval(self.lineEdit5[m].text()))
                m+=1
        m = 0
        for i in range(self.nm):
            for j in range(self.ng):
                self.Vect3[m] = eval(self.lineEdit5[m].text())
                m+=1
        self.close() 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window4(QWidget):
    def __init__(self,ngroup,nmat,order,parent=None):
        super(Window4, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        self.setMinimumSize(QtCore.QSize(100, 100))
        self.ng = ngroup
        self.nm = nmat
        self.ord = order
        lab5  = [0]*ngroup*ngroup*nmat*(order+1)
        self.lineEdit6 = [0]*ngroup*ngroup*nmat*(order+1)
        self.Vect4  = [0]*ngroup*ngroup*nmat*(order+1)
        self.SigS   = [0]*ngroup*ngroup*nmat*(order+1)
        self.tab4   = [[0]*ngroup]*ngroup*nmat*(order+1)
        layout = QGridLayout()
        #------------------------------------
        topLabel5 = QLabel("<font color=blue > Scatter Matrix Cross Section (SigS):</font>") 
        m1 =  0
        m  =  0
        layout.addWidget(topLabel5,0,0)
        for i in range(nmat):
            lab5[m] = QLabel("Material %s" %(i+1))
            layout.addWidget(lab5[m], m1+1, 0)
            for k in range(ngroup):
                lab5[k] = QLabel("G %s" %(k+1))
                layout.addWidget(lab5[k], m1+2, k+2)
            for l in range(order+1):
                lab5[k] = QLabel("Legendre Order L = %s" %(l))
                layout.addWidget(lab5[k], m1+2, 0)
                for j in range(ngroup):
                    lab5[j] = QLabel("G %s" %(j+1))
                    layout.addWidget(lab5[j], m1+3, 1)
                    for n in range(ngroup):
                        self.lineEdit6[m] = QLineEdit()
                        layout.addWidget(self.lineEdit6[m], m1+3, n+2)
                        self.lineEdit6[m].insert(str(self.Vect4[m]))
                        m = m + 1
                    m1 = 1+m1
                m1 = 1+m1
            m1 = 1+m1
        #créer un bouton
        self.bouton = QPushButton(u"Save")
        self.bouton.clicked.connect(self.save4)
        layout.addWidget(self.bouton, m1+1, 0)
        self.setLayout(layout)

    def save4(self):
        self.SigS = []
        m=0
        for i in range(self.nm):
            self.SigS.append([])
            for l in range(self.ord+1):
                self.SigS[i].append([])
                for j in range(self.ng):
                    self.SigS[i][l].append([])
                    for n in range(self.ng):
                        self.SigS[i][l][j].append(eval(self.lineEdit6[m].text()))
                        m+=1

        m = 0
        for i in range(self.nm):
            for l in range(self.ord+1):
                for j in range(self.ng):
                    for n in range(self.ng):
                        self.Vect4[m] = eval(self.lineEdit6[m].text())
                        m +=1
        self.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Window5(QWidget):
    def __init__(self,ngroup,nmat,parent=None):
        super(Window5, self).__init__(parent)
        self.setWindowTitle("Insert Input Parameters")
        self.setWindowFlags(QtCore.Qt.CustomizeWindowHint | QtCore.Qt.WindowCloseButtonHint)
        self.setMinimumSize(QtCore.QSize(100, 100))
        self.ng = ngroup
        self.nm = nmat
        self.Vect5  = [0]*ngroup*nmat
        self.Chi  = [0]*ngroup*nmat
        self.tab5  = [[0]*ngroup]*nmat
        self.lineEdit7 = [0]*ngroup*nmat
        lab4  = [0]*ngroup*nmat
        layout = QGridLayout()
        #--------------------------------------------
        topLabel1 = QLabel("<font color=blue > Density Function for Neutrons (Chi):</font>") 
        m = 0
        layout.addWidget(topLabel1,0,0)
        for j in range(nmat):
            lab4[j] = QLabel("Material %s" %(j+1))
            layout.addWidget(lab4[j], j+2, 0)
            for i in range(ngroup):
                lab4[i] = QLabel("G %s" %(i+1))
                layout.addWidget(lab4[i], 1, i+1)
                self.lineEdit7[m] = QLineEdit()
                layout.addWidget(self.lineEdit7[m], j+2, i+1)
                self.lineEdit7[m].insert(str(self.Vect5[m]))
                m=m+1
        #réer un Bouton
        self.bouton = QPushButton(u"Save")
        self.bouton.clicked.connect(self.save5)
        layout.addWidget(self.bouton, j+3, 0)
        self.setLayout(layout)

    def save5(self):
        m=0
        self.Vect5 = [0]*self.nm*self.ng
        del self.Chi[:]
        for i in range(self.nm):
            self.Chi.append([])
            for j in range(self.ng):
                self.Chi[i].append(eval(self.lineEdit7[m].text()))
                m+=1
        m = 0
        for i in range(self.nm):
            for j in range(self.ng):
                self.Vect5[m] = eval(self.lineEdit7[m].text())
                m+=1
        self.close() 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# An Example application QWidget containing the textEdit_4 to redirect stdout to
class Application(QtWidgets.QMainWindow):

    def __init__(self, title= "Default", parent=None):
        super(Application, self).__init__(parent)
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
        self.ui.pushButton_2.clicked.connect(self.visualization)
        self.ui.pushButton_12.clicked.connect(self.visualization2d)
        self.ui.pushButton_3.clicked.connect(self.plot)
        self.ui.pushButton_4.clicked.connect(self.compile)
        self.ui.pushButton_18.clicked.connect(self.compile)
        self.ui.pushButton_5.clicked.connect(self.run)
        self.ui.pushButton_19.clicked.connect(self.powerpf)
        self.ui.pushButton_6.clicked.connect(self.power)
        self.ui.pushButton_20.clicked.connect(self.run)
        self.ui.pushButton.clicked.connect(self.new1)
        self.ui.pushButton_7.clicked.connect(self.new2)
        self.ui.pushButton_8.clicked.connect(self.new3)
        self.ui.pushButton_9.clicked.connect(self.new4)
        self.ui.pushButton_26.clicked.connect(self.new5)
        self.ui.pushButton_23.clicked.connect(self.new6) # chi (14) nusigf(15) sigt(17) sigs(13) size(23) data (16)
        #self.ui.pushButton_17.clicked.connect(self.new7)
        #self.ui.pushButton_15.clicked.connect(self.new8)
        #self.ui.pushButton_13.clicked.connect(self.new9)
        #self.ui.pushButton_14.clicked.connect(self.new10)
        self.ui.pushButton_11.clicked.connect(self.data_up)
        self.ui.pushButton_16.clicked.connect(self.data_up2D)
        self.statusbar = self.statusBar()
        self.ui.textEdit_4.cursorPositionChanged.connect(self.cursorPosition_2)
        self.ui.textEdit_3.cursorPositionChanged.connect(self.cursorPosition)
        self.ui.radioButton_4.toggled.connect(lambda:self.methods_1D(self.ui.radioButton_4))
        self.ui.radioButton_2.toggled.connect(lambda:self.methods_1D(self.ui.radioButton_2))
        self.ui.radioButton_3.toggled.connect(lambda:self.methods_1D(self.ui.radioButton_3))  
        self.ui.radioButton.toggled.connect(lambda:self.methods_2D(self.ui.radioButton)) 

        sys.stdout = EmittingStream(textWritten=self.normalOutputWritten)
        sys.stderr = EmittingStream(textWritten=self.normalOutputWritten)
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
        
        M00 = open('app/link/script00.py', "r" ).read() 
        if M00 == 'CP1D':
            self.ui.radioButton_2.setChecked(True) 
        elif M00 == 'SN1D':
            self.ui.radioButton_3.setChecked(True)
        elif M00 == 'MOC1D':
            self.ui.radioButton_4.setChecked(True)       
        elif M00 == 'SN2D':
            self.ui.radioButton.setChecked(True) 

        # Multi-group Cross Sections
        self.ui.comboBox_12.currentIndexChanged.connect(self.Cross)
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
            self.ui.comboBox_12.setCurrentText('Chi')

        # Constructive geometry
        self.ui.comboBox_7.currentIndexChanged.connect(self.ConstGeom)
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def new1(self):
        nregion  = self.ui.spinBox_2.value()
        nmat     = self.ui.spinBox_3.value()  
        if  nregion == 0:
            QMessageBox.warning(self, "Warning", "Enter the Total Number of Regions")
        elif nmat == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of materials")
        elif nmat > nregion:
            QMessageBox.warning(self, "Warning", "the number of materials must not exceed the number of regions")
        else:
            self.wind1 = Window1(nregion,nmat)
            self.wind1.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def new2(self):
        ngroup   = self.ui.spinBox.value()
        nmat     = self.ui.spinBox_3.value() 
        if  ngroup == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of energy group")
        elif nmat == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of materials")
        else:
            self.wind2 = Window2(ngroup,nmat)
            self.wind2.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def new3(self):
        ngroup   = self.ui.spinBox.value()
        nmat     = self.ui.spinBox_3.value()  
        if  ngroup == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of energy group")
        elif nmat == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of materials")
        else:
            self.wind3 = Window3(ngroup,nmat)
            self.wind3.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def new4(self):
        ngroup   = self.ui.spinBox.value()
        nmat     = self.ui.spinBox_3.value()  
        order    = self.ui.spinBox_5.value() 
        if  ngroup == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of energy group")
        elif nmat == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of materials")
        else:
            self.wind4 = Window4(ngroup,nmat,order)
            self.wind4.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# self.ui.spinBox.value() Energy Group Number 1D
# self.ui.spinBox_8.value() Energy Group Number 2D
# self.ui.spinBox_7.value() Tolerance 1D
# self.ui.spinBox_6.value() Maximum Iteration 1D
# self.ui.spinBox_2.value() Number of Regions 1D
# self.ui.spinBox_3.value()  Number of Materials 1D
# self.ui.spinBox_4.value() Angular Discretizations 1D
# self.ui.spinBox_5.value() Legendre Order 1D
# self.ui.spinBox_9.value() Number of Pin Cells 2D
# self.ui.spinBox_10.value() Number of Materials 2D
# self.ui.spinBox_11.value() Number of Assemblies
# self.ui.spinBox_16.value() Tolerance 2D
# self.ui.spinBox_15.value() Maximum Iteration 2D
# self.ui.spinBox_12.value() Angular Discretizations 2D
# self.ui.spinBox_13.value() Legendre Order 2D
# self.ui.spinBox_14.value() Pin x-y Number
# Function for Multi-group Cross Sections
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
            elif M06 == 'Chi':
                self.wind13 = Window13(ngroup,nmat)
                self.wind13.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def new7(self):
        ngroup   = self.ui.spinBox_8.value()
        nmat     = self.ui.spinBox_10.value() 
        if  ngroup == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of energy group")
        elif nmat == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of materials")
        else: 
            self.wind2 = Window2(ngroup,nmat)
            self.wind2.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def new8(self):
        ngroup   = self.ui.spinBox_8.value()
        nmat     = self.ui.spinBox_10.value()  
        if  ngroup == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of energy group")
        elif nmat == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of materials")
        else:
            self.wind3 = Window3(ngroup,nmat)
            self.wind3.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def new9(self):
        ngroup   = self.ui.spinBox_8.value()
        nmat     = self.ui.spinBox_10.value()  
        order    = self.ui.spinBox_13.value() 
        if  ngroup == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of energy group")
        elif nmat == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of materials")
        else:
            self.wind4 = Window4(ngroup,nmat,order)
            self.wind4.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def new10(self):
        ngroup   = self.ui.spinBox_8.value()
        nmat     = self.ui.spinBox_10.value()  
        if  ngroup == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of energy group")
        elif nmat == 0:
            QMessageBox.warning(self, "Warning", "Enter the number of materials")
        else:
            self.wind5 = Window5(ngroup,nmat)
            self.wind5.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def data_up(self):
        Geometry_type = open('app/link/script01.py', "r" ).read()
        ngroup  = self.ui.spinBox.value()
        nregion = self.ui.spinBox_2.value()
        nmat    = self.ui.spinBox_3.value()
        ngauss  = self.ui.spinBox_4.value()
        order   = self.ui.spinBox_5.value()
        max_it  = self.ui.spinBox_6.value()
        Tol     = self.ui.spinBox_7.value()
        if  ngroup == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres Energy Group Number")
        elif nregion == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres Number of Regions")
        elif nmat == 0:
            QMessageBox.warning(self, "Warning", "Enter the Input Parametres Number of Materials")
        else: 
            wind1 = Window1(nregion,nmat) 
            wind2 = Window2(ngroup,nmat)
            wind3 = Window3(ngroup,nmat)
            wind4 = Window4(ngroup,nmat,order)
            wind5 = Window5(ngroup,nmat)
            filename = open("app/input/input.json",'w')
            open('app/link/script.py', "w" ).write(os.path.abspath(os.path.dirname( __file__)) +'/input/input.json')
            filename.write('{ \n  "data": { \n    "parameter": { \n      "id": 100,')
            filename.write('\n      "Total number of energy groups": '+ str(ngroup) +',')
            filename.write('\n      "Total number of Materials": '+ str(nmat) +',')
            filename.write('\n      "Total number of regions": '+ str(nregion) +',')
            filename.write('\n      "Which material fills each region": '+ str(self.wind1.REGMAT) +',')
            if Geometry_type in ['Cylindrical Geometry','Spherical Geometry']:
                filename.write('\n      "Ray for each region per [cm]": '+ str(self.wind1.Delta) +',')
            else:
                filename.write('\n      "Size of each region [cm]": '+ str(self.wind1.Delta) +',')
            filename.write('\n      "Number of fine meshes": '+ str(self.wind1.NFMR) +',')  
            filename.write('\n      "Number of Angular Discretization": '+ str(ngauss) +',')
            filename.write('\n      "The l-order Legendre polynomial": '+ str(order) +',')
            filename.write('\n      "Maximum Number of Iterations": '+ str(max_it) +',')
            filename.write('\n      "Criterion of Keff convergence": '+ '1.0E'+str(Tol))
            filename.write('\n    }, \n    "materials": [')
            # Ici Boucle
            for i in range(nmat):
                filename.write('\n      { \n        "id": '+ str(i+1) +', \n        "nom": "material '+ str(i+1) +'",') 
                filename.write('\n        "XSTotal": ' + str(self.wind2.SigT[i][:]) +','
                                   + '\n        "XSNuFission": '+ str(self.wind3.NuSigF[i][:])+',')
                filename.write('\n        "XSScatter Matrix":'+str(self.wind4.SigS[i][:][:][:])+','+ 
                 '\n        "XSChi":  '+str(self.wind5.Chi[i][:]))
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
    def data_up2D(self):   
# self.ui.spinBox.value() Energy Group Number 1D
# self.ui.spinBox_8.value() Energy Group Number 2D
# self.ui.spinBox_7.value() Tolerance 1D
# self.ui.spinBox_6.value() Maximum Iteration 1D
# self.ui.spinBox_2.value() Number of Regions 1D
# self.ui.spinBox_3.value()  Number of Materials 1D
# self.ui.spinBox_4.value() Angular Discretizations 1D
# self.ui.spinBox_5.value() Legendre Order 1D
# self.ui.spinBox_9.value() Number of Pin Cells 2D
# self.ui.spinBox_10.value() Number of Materials 2D
# self.ui.spinBox_11.value() Number of Assemblies
# self.ui.spinBox_16.value() Tolerance 2D
# self.ui.spinBox_15.value() Maximum Iteration 2D
# self.ui.spinBox_13.value() Angular Discretizations 2D
# self.ui.spinBox_12.value() Legendre Order 2D
# self.ui.spinBox_14.value() Pin x-y Number
        ngroup   = self.ui.spinBox_8.value()
        np       = self.ui.spinBox_9.value()   
        nmat     = self.ui.spinBox_10.value() 
        na       = self.ui.spinBox_11.value() 
        nxy      = self.ui.spinBox_14.value()
        nxya     = self.ui.spinBox_17.value()
        ordinate = self.ui.spinBox_12.value()
        order    = self.ui.spinBox_13.value() 
        max_it   = self.ui.spinBox_15.value()
        Tol      = self.ui.spinBox_16.value()
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
        else: 
            #wind6 = Window1(nx,ny) 
            #wind2 = Window2(ngroup,nmat)
            #wind3 = Window3(ngroup,nmat)
            #wind4 = Window4(ngroup,nmat,order)
            #wind5 = Window5(ngroup,nmat)
            wind6  = Window6(nxy,nxy,np)
            wind7  = Window7(nxya,na)
            wind8  = Window8(na)
            wind9  = Window9(ngroup,nmat) 
            wind10 = Window10(ngroup,nmat)
            wind11 = Window11(ngroup,nmat)
            wind12 = Window12(ngroup,nmat)
            wind13 = Window13(ngroup,nmat)
            filename = open("app/input/input.json",'w')
            open('app/link/script.py', "w" ).write(os.path.abspath(os.path.dirname( __file__)) +'/input/input.json')
            filename.write('{ \n  "data": { \n    "parameter": { \n      "id": 100,')
            filename.write('\n      "Total number of energy groups": '+ str(ngroup) +',')
            filename.write('\n      "Total number of materials": '+ str(nmat) +',')
            filename.write('\n      "Total number of pin cells": '+ str(np) +',')
            filename.write('\n      "Total number of assemblies": '+ str(na) +',')
            filename.write('\n      "Core": '+ str(self.wind8.core) +',')
            filename.write('\n      "Number of angular discretizations": '+ str(ordinate) +',')
            filename.write('\n      "The l-order Legendre polynomial": '+ str(order) +',')
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
    @QtCore.pyqtSlot()
    def WorkReply(self):
        """Exécuté lorsque le processus envoie des infos à afficher.
           La chaine renvoyée par data() est de type byte, terminée
           par une fin de ligne. 
           L'encodage dépend de la commande lancée.
        """
        data = self.process.readAllStandardOutput().data()
        ch = str(data, encoding="utf-8").rstrip()
        print(ch)
 
    @QtCore.pyqtSlot()
    def WorkFinished(self):
        """exécuté à la fin du processus
        """
        if self.process!=None:
            # le processus vient de se terminer: on fait le ménage
            self.process.readyReadStandardOutput.disconnect()
            self.process.finished.disconnect()
            QMessageBox.information(self, "Message", "Compilation process is completed.")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    @QtCore.pyqtSlot()
    def run(self):
        """exécuté lors du clic sur le bouton
        """
        self.process = QtCore.QProcess()
        self.process.setProcessChannelMode(QtCore.QProcess.MergedChannels)
        self.process.readyReadStandardOutput.connect(self.WorkReply)
        self.process.finished.connect(self.WorkFinished)
        commande = "python run.py"
        self.process.start(commande, QtCore.QIODevice.ReadWrite)
        self.process.waitForStarted()
 
    @QtCore.pyqtSlot()
    def WorkReply(self):
        """Exécuté lorsque le processus envoie des infos à afficher.
           La chaine renvoyée par data() est de type byte, terminée
           par une fin de ligne. 
           L'encodage dépend de la commande lancée.
        """
        data = self.process.readAllStandardOutput().data()
        ch = str(data).rstrip()
        print(ch)
 
    @QtCore.pyqtSlot()
    def WorkFinished(self):
        """exécuté à la fin du processus
        """
        if self.process!=None:
            # le processus vient de se terminer: on fait le ménage
            self.process.readyReadStandardOutput.disconnect()
            self.process.finished.disconnect()
            #QMessageBox.warning(self, "Warning", "Running case finished.")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def visualization(self):
        Geometry_type = open('app/link/script01.py', "r" ).read()
        if Geometry_type in ['Slab Geometry']:
            filename = open('app/link/script.py', "r" ).read()
            with open(filename) as json_data:
                data = json.load(json_data)
                l = data['data']['parameter']['Size of each region [cm]']
                self.nmat = data['data']['parameter']['Total number of Materials']
                self.nregion = data['data']['parameter']['Total number of regions'] 
                self.regmat = data['data']['parameter']['Which material fills each region']
            fig, ax = plt.subplots()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
            width = [0]
            som = [0]
            for i in range(self.nregion):
              width.append(l[i])
              som.append(sum(width))
            rectangles = []
            red_patch = []
            m = som[0]

            for n in range(self.nregion):
                rectangles.append(mpatch.Rectangle((som[n],-5), abs(som[n+1])-abs(som[n]), 10))
            colr = ['black','gray','darkgray','silver','gainsboro','whitesmoke','snow','white']
            i=0
            for r in rectangles:
                ax.add_artist(r) 
                r.set_facecolor(color=colr[self.regmat[i]-1])
                rx, ry = r.get_xy()
                cx = rx + r.get_width()/2.0
                cy = ry + r.get_height()/2.0
                i=i+1
            for i in range(self.nmat):
                red_patch.append(mpatch.Patch(color=colr[i], label="Mat %s" %(i+1)))

            ax.set_ylim((-10, 10))
            ax.set_xlim((min(som), max(som)))
            ax.set_xlabel('X [cm]')
            ax.set_title('Color by Materials') 
            plt.legend(handles=red_patch, loc='center left', bbox_to_anchor=(1, 0.5))
            plt.show()
        else:
            filename = open('app/link/script.py', "r" ).read()
            with open(filename) as json_data:
                data = json.load(json_data)
                l = data['data']['parameter']['Ray for each region per [cm]']
                self.nregion = data['data']['parameter']['Total number of regions']
                self.nmat = data['data']['parameter']['Total number of Materials']
                rayon = data['data']['parameter']['Ray for each region per [cm]']
            fig, ax = plt.subplots()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
            red_patch = []
            levels = [0.0]+rayon 
            xlist = np.linspace(-max(l)-2, max(l)+2, 100)
            ylist = np.linspace(-max(l)-2, max(l)+2, 100)
            X, Y = np.meshgrid(xlist, ylist)
            Z = np.sqrt(X ** 2 + Y ** 2 )
            contour = plt.contour(X, Y, Z, levels, colors='k')
            plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=12)
            plt.axes().set_aspect("equal")
            colr = ['r','yellow',"c","Grey",'b', 'y','c', 'm', 'y', 'k', 'w']
            c = ('r','yellow',"c","Grey",'b', 'y','c', 'm', 'y', 'k', 'w')
            contour_filled = plt.contourf(X, Y, Z, levels, colors=c)
            for i in range(self.nmat):
                red_patch.append(mpatch.Patch(color=colr[i], label="Mat %s" %(i+1)))
            plt.title('Color by Materials')
            plt.xlabel('R [cm]')
            plt.ylabel('R [cm]')
            plt.legend(handles=red_patch, loc='center left', bbox_to_anchor=(1, 0.5))
            plt.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def visualization2d(self):
        filename = open('app/link/script.py', "r" ).read()
        fmm_id =[[]]
        assembly = []
        pin = []
        regmat = []
        nom    = []
        width_x = []
        width_y = []
        with open(filename) as json_data:
            data = json.load(json_data)
            nmat = data['data']['parameter']['Total number of materials']
            npc = data['data']['parameter']['Total number of pin cells']  # Number of Pin Cell
            na  = data['data']['parameter']['Total number of assemblies'] 
            core = data['data']['parameter']['Core']

            for i in range(na):
                assembly.append(data['data']['Assemblies'][i]['assembly'])

            for i in range(npc):
                pin.append(data['data']['PinCells'][i]['mat_fill'])
            Height = len(core[0])*len(pin[0][0])*len(assembly[0])
            width = len(core[0])*len(pin[0][0])*len(assembly[0])
            fmm_id =numpy.zeros(((Height),(width)),dtype='i')

            i1=0
            for j in range(len(core[0])):
                for k in range(len(assembly[0])):
                    for m in range(len(pin[0])):
                        i2=0
                        for i in range(len(core[0])):
                            for l in range(len(assembly[0])):
                                for n in range(len(pin[0])):
                                    fmm_id[i1][i2] = pin[ assembly[ core[j][i] -1][k][l] -1][m][n]
                                    i2+=1
                        i1+=1
            

            nx = len(np.amax(data['data']['parameter']['Core'],axis=0))
            ny = len(np.amax(data['data']['parameter']['Core'],axis=1))
            
            nxx = len(assembly[0])
            nyy = len(assembly[0])
            NX = nxx*nx
            NY = ny*nyy
            width_x.append(data['data']['PinCells'][0]['width_x'])
            width_y.append(data['data']['PinCells'][0]['width_y'])
            nx =  len(width_x[0])*NX
            ny =  len(width_y[0])*NY
            xcm  = width_x[0]*NX
            ycm  = width_y[0]*NY
            for j in range(nmat):
                nom.append(data['data']['materials'][j]['nom'])
            for i in range(npc):
                regmat.append(data['data']['PinCells'][i]['mat_fill'])

            fig, ax = plt.subplots()
            widthx = [0]
            widthy = [0]
            somx = [0]
            somy = [0]
            for i in range(nx):
                widthx.append(xcm[i])
                somx.append(sum(widthx))
            for j in range(ny):
                widthy.append(ycm[j])
                somy.append(sum(widthy))
            rectangles = []
            red_patch = []
            mx = somx[0]
            my = somy[0]

            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.77, box.height])
            for i in range(nx):
                for j in range(ny):
                    rectangles.append(mpatch.Rectangle((somx[i],somy[j]), xcm[i], ycm[j],linewidth=0.0,edgecolor='k'))
            colr = ['black','gray','darkgray','silver','gainsboro','whitesmoke','snow','white']
            #colr = ['black','dimgrey','grey','darkgrey','lightgrey','whitesmoke','white','b']
            for i in range(nmat):
                red_patch.append(mpatches.Patch(color=colr[i], label=nom[i]))
            n=nx-1
            m=0
            for r in rectangles:
                ax.add_artist(r) 
                r.set_facecolor(color=colr[int(fmm_id[n][m])-1])

                rx, ry = r.get_xy()
                cx = rx + r.get_width()/2.0
                cy = ry + r.get_height()/2.0
                n-=1
                if (n<0):
                    n=nx-1
                    m+=1

            ax.set_ylim((min(somy), max(somy)))
            ax.set_xlim((min(somx), max(somx)))
            #ax.set_xticklabels([])
            #ax.set_yticklabels([]) 
            ax.set_xlabel('X [cm]')
            ax.set_ylabel('Y [cm]')
            #ax.set_title('Color by Materials') 
            clb = plt.legend(handles=red_patch, loc='center left', title="Materials", fontsize='small', bbox_to_anchor=(1, 0.5))
            plt.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def P_Ordinate(self):
        file_name = 'app/Output/Ordinates.h'
        positions = np.loadtxt(file_name)
        ss = len(positions)
        red_patch = []
        if (ss == 4):
            x = np.array( [positions[0]] )
            y = np.array( [positions[1]] )
            z = np.array( [positions[2]] )
            w = np.array( [positions[3]] )
        else:
            x = positions[:,0]
            y = positions[:,1]
            z = positions[:,2]
            w = positions[:,3]
        ww = []
        i=1
        for a in w:
            if a not in ww:
                ww.append(a)
                red_patch.append(mpatch.Patch(label="W%s= %s" %(i,a)))
                i+=1

        w = 50. * w / w.max()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim3d([0, 1])
        ax.set_ylim3d([0, 1])
        ax.set_zlim3d([0, 1])
        ax.set_autoscale_on(False)
        n = 100
        zr = []
        for i in range(0, len(x)):
            add = True
            for j in range(0, len(zr)):
                if (zr[j][0] == z[i]):
                    add = False
            if (add==True):
                zr.append( [z[i], ((x[i]**2.)+(y[i]**2.))**0.5] )
 
        circle = []
        #box = Bbox(((-1,-1),(1,1)))

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        for i in range(0, len(zr)):
            r = zr[i][1]
            d = 2*r
            zed = zr[i][0]
            circle.append(Arc((0,0), width=d, height=d, angle=0.0, theta1=0.0, theta2=90.0, color="white", alpha=0.2))
            ax.add_patch(circle[-1])
            art3d.pathpatch_2d_to_3d(circle[-1], z=zed, zdir="x")
            circle[-1].set_clip_on(True)
    
            circle.append(Arc((0,0), width=d, height=d, angle=0.0, theta1=0.0, theta2=90.0, color="white", alpha=0.2))
            ax.add_patch(circle[-1])
            art3d.pathpatch_2d_to_3d(circle[-1], z=zed, zdir="y")
            circle[-1].set_clip_on(True)
    
            circle.append(Arc((0,0), width=d, height=d, angle=0.0, theta1=0.0, theta2=90.0, color="white", alpha=0.2))
            ax.add_patch(circle[-1])
            art3d.pathpatch_2d_to_3d(circle[-1], z=zed, zdir="z")
            circle[-1].set_clip_on(True)
        for i in range(0, len(x)):
            custom = ( x[i], y[i], z[i] )
            ax.scatter(x[i],y[i],z[i],color=custom,s=w[i]*5,zorder=2)
        ax.view_init(elev=20., azim=25)
        ax.set_xlabel(u'\u00B5')
        ax.set_ylabel(u'\u03B7') 
        ax.set_zlabel(u'\u03BE')
        if (len(positions) == 4):
            tt = '2'
        else:
            delta = 4 +4*8*len(positions)
            tt = str(int((-1+ math.sqrt(delta))/2))
        plt.title(r'$S_{'+''.join([i for i in tt[0:]])+'}$')
        #fig.savefig(file_name+".png", dpi=100, facecolor='white')
        plt.legend(handles=red_patch,handlelength=False, handletextpad=0, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def plot(self):
        filename = open('app/link/script.py', "r" ).read()
        with open(filename) as json_data:
            data = json.load(json_data)
            ng = data['data']['parameter']['Total number of energy groups']
         
        data = []
        M00 = open('app/link/script00.py', "r" ).read()
        if M00 == 'CP1D':
            data = np.loadtxt('app/Output/flux_cp.h')
        elif M00 == 'SN1D':
            data = np.loadtxt('app/Output/flux_sn.h')
            #N = data['data']['parameter']['Number of Angular Discretization']
        elif M00 == 'MOC1D':
            data = np.loadtxt('app/Output/flux_moc.h')
        elif M00 == 'SN2D':
            if  ng == 0:
                QMessageBox.warning(self, "Warning", "Enter the Input Parametres Energy Group Number")
                return
            else:
                data = np.loadtxt('app/Output/FG'+str(ng))
        else:
            QMessageBox.warning(self, "Warning", "select the calculation method.")
        if int(len(data)) >= 0:           
            if M00 != 'SN2D':
                fig, ax = plt.subplots()
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
                matrix = []  
                for line in data:
                    matrix.append(line)
                max_columns = len(matrix[0]) - 1
                max_rows = len(matrix)
                x = [matrix[rownum][0] for rownum in range(max_rows)]
                y = [[matrix[rownum][colnum + 1] for rownum in range(max_rows)] for colnum in range(max_columns)]
                p = [0]*max_columns
                for i in range(max_columns):
                    key = 'Group', int(i+1)
                    p[i] = plt.plot(x,y[i] ,label="Group %s" %(max_columns-i),linewidth=1)
                if M00 == 'CP1D':
                    plt.title("CP Method")
                elif M00 == 'SN1D':
                    plt.title("SN Method") 
                elif M00 == 'MOC1D':
                    plt.title("MOC")  
                plt.xlabel('Distance [cm]')
                plt.ylabel('Normalized Flux')
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                plt.show()
            elif M00 == 'SN2D':        
                max_columns = len(data[0]) - 1
                max_rows = len(data)
                x = [data[rownum+1][0] for rownum in range(max_rows-1)] 
                y = [data[0][colnum + 1] for colnum in range(max_columns)]
                z = [[data[rownum+1][colnum + 1] for rownum in range(max_rows-1)] for colnum in range(max_columns)]
                x = np.array(x)
                y = np.array(y)
                z = np.array(z)
                x,y = np.meshgrid(x, y, sparse=True)
                fig = matplotlib.pyplot.figure()

                ax = fig.add_subplot(111)
                ax.set_xlabel('X [cm]')                         
                ax.set_ylabel('Y [cm]') 
                ax.set_title('Energy Group '+str(ng))
                #ax.set_title(r'$S_{'+''.join(str(N))+'}$, Energy Gr '+str(ng))
                im = ax.imshow(z, interpolation='bilinear', cmap='jet', 
                origin='lower', extent=[0, abs(x).max(), 0,abs(y).max()])
                clb = fig.colorbar(mappable=im, format='%.1E')
                clb.set_label('Normalized scalar flux')
                clb.update_normal(im) 


                callback = Index(ng,ax,clb)
                axprev = plt.axes([0.785, 0.01, 0.05, 0.06])
                axnext = plt.axes([0.85, 0.01, 0.05, 0.06])
                bnext = Button(axnext, '>')
                bnext.on_clicked(callback.next)
                bprev = Button(axprev, '<')
                bprev.on_clicked(callback.prev)
                plt.show()
        else:
            QMessageBox.warning(self, "Warning", "Select More than a Fine Number of Meshes")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def powerpf(self):
        data = np.loadtxt('app/Output/PowerPeakingFactor.h')
        import matplotlib.pyplot
        max_columns = len(data[0]) - 1
        max_rows = len(data)
        x = [data[rownum+1][0] for rownum in range(max_rows-1)] 
        y = [data[0][colnum + 1] for colnum in range(max_columns)]
        z = [[data[rownum+1][colnum + 1] for rownum in range(max_rows-1)] for colnum in range(max_columns)]
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        x,y = np.meshgrid(x, y)

        fig = matplotlib.pyplot.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('X [cm]')                         
        ax.set_ylabel('Y [cm]') 
        im = ax.imshow(z, interpolation='bilinear', cmap='jet', 
             origin='lower', extent=[0, abs(x).max(), 0,abs(y).max()])
        clb=fig.colorbar(im, format='%.1E')
        clb.set_label('Normalized Pin Power Distribution')
        plt.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def power(self):
        data = np.loadtxt('app/Output/PowerDistribution.h')
        max_columns = len(data[0]) - 1
        max_rows = len(data)
        x = [data[rownum+1][0] for rownum in range(max_rows-1)] 
        y = [data[0][colnum + 1] for colnum in range(max_columns)]
        z = [[data[rownum+1][colnum + 1] for rownum in range(max_rows-1)] for colnum in range(max_columns)]
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        x,y = np.meshgrid(x, y)

        fig = matplotlib.pyplot.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('X [cm]')                         
        ax.set_ylabel('Y [cm]') 
        im = ax.imshow(z, interpolation='bilinear', cmap='jet', 
             origin='lower', extent=[0, abs(x).max(), 0,abs(y).max()])
        clb = fig.colorbar(im, format='%.1E')
        #ax.set_xticklabels([])
        #ax.set_yticklabels([])
        #clb.ax.set_title('Peaking Factor Distribution ', loc='left', fontsize=9)
        #plt.xticks(rotation=45)
        #plt.yticks(rotation=45)
        #ax.xaxis.tick_top()
        clb.set_label('Power Distribution')
        plt.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def aboutQt(self):
        self.infoLabel.setText("Invoked <b>Help|About Qt</b>") 
    def Cross(self,b):
        open('app/link/script11.py', "w" ).write(self.ui.comboBox_12.currentText()) 
    def ConstGeom(self,b):
        open('app/link/script10.py', "w" ).write(self.ui.comboBox_7.currentText()) 
    def geometry(self,b):
        open('app/link/script01.py', "w" ).write(self.ui.comboBox.currentText()) 	
        #print (b,self.ui.comboBox.currentText())
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
        open('app/link/script.py', "w" ).write(str(self.filename))
        if self.filename != "":
            self.isFileOpen = True
            file = open(self.filename,"r")
            self.setWindowTitle(self.title + ':' + self.filename)
            self.ui.textEdit_4.show()
            self.ui.textEdit_4.setText(file.read())
            file.close()

    def saveFile(self):
        import time
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
        import time
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

    #@pyqtSlot(str)
    #def append_text(self,text):
        #self.ui.textEdit_3.moveCursor(QTextCursor.End)
        #self.ui.textEdit_3.insertPlainText( text )


    @pyqtSlot()
    def start_thread1(self):
        self.thread = QThread()
        self.long_running_thing = LongRunningThing()
        self.long_running_thing.moveToThread(self.thread)
        self.thread.started.connect(self.long_running_thing.run)
        self.thread.start()
        self.thread.quit()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    @QtCore.pyqtSlot()
    def compile(self):
        """exécuté lors du clic sur le bouton
        """
        #print ('Connecting process')
        self.process = QtCore.QProcess(self)
        self.process.readyReadStandardOutput.connect(self.stdoutReady)
        self.process.readyReadStandardError.connect(self.stderrReady)
        #print ('Starting process')
        self.process.start('python', ['compile.py'])



    @pyqtSlot()
    def run(self):
        #print ('Connecting process')
        self.process = QtCore.QProcess(self)
        self.process.readyReadStandardOutput.connect(self.stdoutReady)
        self.process.readyReadStandardError.connect(self.stderrReady)
        #print ('Starting process')
        self.process.start('python', ['app/python.py'])
        

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





# Create Queue and redirect sys.stdout to this queue
#queue = Queue()
#sys.stdout = WriteStream(queue)

#redirect stdout

qapp = QApplication(sys.argv)  
app = Application(u'OpenNTP')
app.show()
qapp.exec_()







