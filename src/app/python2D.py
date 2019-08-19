#! /usr/bin/env python
#! -*- coding:utf-8 -*-
import SlabSN2D
from matplotlib import pyplot as plt
import numpy as np
import json
from datetime import datetime
import subprocess
start = datetime.now()
sigtt  = []
sigss  = []
sigt   = []
nusigf = []
sigs   = []
chi    = []
vol    = []
e = "                         _   _ _____ ____       _____ ____  ____  _   _ "
a = "   ___  _ __   ___ _ __ | \ | |_   _|  _ \     | ____|  _ \/ ___|| \ | |"
b = "  / _ \| '_ \ / _ | '_ \|  \| | | | | |_) _____|  _| | |_) \___ \|  \| |"
c = " | (_) | |_) |  __| | | | |\  | | | |  __|_____| |___|  _ < ___) | |\  |" 
d = "  \___/| .__/ \___|_| |_|_| \_| |_| |_|        |_____|_| \_|____/|_| \_|"
f = "       |_|                                                              "
print e
print a
print b
print c
print d
print f
filename      = open('app/link/script.py', "r" ).read()
Geometry_type = open('app/link/script01.py', "r" ).read()
BC            = open('app/link/script02.py', "r" ).read()
scheme1       = open('app/link/script03.py', "r" ).read()
scheme2       = open('app/link/script04.py', "r" ).read()
Methods_2D    = open('app/link/script00.py', "r" ).read()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
if Methods_2D == 'SN':
#=========================================================================================
#  Loading data for SN Method 2D
#========================================================================================= 
    regmat    = []
    nfmesh_xy = []
    with open(filename) as json_data:
        data = json.load(json_data)
        ng = data['data']['parameter']['Total number of energy groups'] 
        nmat = data['data']['parameter']['Total number of Materials']
        nreg_x = data['data']['parameter']['Total number of X regions']
        nreg_y = data['data']['parameter']['Total number of Y regions'] 
        cell_x = data['data']['parameter']['X region thickness per [cm]']
        cell_y = data['data']['parameter']['Y region thickness per [cm]']
        for i in range(nreg_y):
            regmat.append(data['data']['parameter']['Which material goes in each cell'][i])
            nfmesh_xy.append(data['data']['parameter']['XY number of fine meshes in each cell'][i])
        ngauss = data['data']['parameter']['Number of Angular Discretization']
        Max_it = data['data']['parameter']['Maximum Iteration']
        order =  data['data']['parameter']['The l-order Legendre polonomial']
        order = order + 1
        eps = data['data']['parameter']['Epsilon Keff']
        for i in range(nmat):
            sigtt.append(data['data']['materials'][i]['XSTotal'])
            nusigf.append(data['data']['materials'][i]['XSNuFission'])
            sigss.append(data['data']['materials'][i]['XSScatter Matrix'])
            chi.append(data['data']['materials'][i]['XSChi']) 
        totnfm_x = sum(np.amax(nfmesh_xy,axis=0)) # x
        totnfm_y = sum(np.amax(nfmesh_xy,axis=1)) # y
        print 'Test',ngauss
        mup,etap,psip,pwi = SlabSN2D.quadrature_set(ngauss)
        print 'coco',mup
        #SlabN2D.fmm_idd(nreg_x)
        #fmmid_xy = SlabN2D.fmm_id(nfmesh_xy,regmat,totnfm_y,totnfm_x,[nreg_x,nreg_y])
        
        #print fmmid_xy
        dx,dy = SlabSN2D.delta_f(totnfm_x,totnfm_y,nfmesh_xy,cell_x,cell_y,[nreg_x,nreg_y])
        mup,etap,psip,pwi = SlabSN2D.quadrature_set(ngauss)
        print mup
        #print mup,etap,psip,pwi 
        #b = SlabSN2D.matrix_axayb(ngauss,mup,etap,psip,fmmid_xy,sigtt,dx,
        #                          dy,[ng,nmat,totnfm_x,totnfm_y,nreg_y,nreg_x])
        #print b
        
        #SlabSN2D.quadrature_set(ngauss)
        #a = SlabSN2D.plot_plm(order-1,order-1)

