#! /usr/bin/python3
#! -*- coding:utf-8 -*-
import app.SlabMOC
import app.SlabSN
import app.CP1D
import app.SN2D
import numpy as np
import json
from datetime import datetime
import subprocess
import os
start  = datetime.now()
sigtt  = []
sigss  = []
sigt   = []
nusigf = []
sigf   = []
sigs   = []
chi    = []
vol    = []
BC2D   = []
regmat  = []
nfmesh  = []
width   = []
ray   = []
assembly  = []
core      = []
albedo = []
PATH          = open(os.getcwd()+'/app/link/script.dir', "r" ).read()
Methods       = open(os.getcwd()+'/app/link/script00.py', "r" ).read()
Geometry_type = open(os.getcwd()+'/app/link/script01.py', "r" ).read()
BC            = open(os.getcwd()+'/app/link/script02.py', "r" ).read()
scheme1       = open(os.getcwd()+'/app/link/script03.py', "r" ).read()
scheme2       = open(os.getcwd()+'/app/link/script04.py', "r" ).read()
bc1           = open(os.getcwd()+'/app/link/script06.py', "r" ).read()
bc2           = open(os.getcwd()+'/app/link/script07.py', "r" ).read()
bc3           = open(os.getcwd()+'/app/link/script08.py', "r" ).read()
bc4           = open(os.getcwd()+'/app/link/script09.py', "r" ).read()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
#==========================================================================================
#  Loading data for MOC 1D
# ========================================================================================= 
if Methods == 'MOC1D':
    with open(PATH) as json_data:
        data = json.load(json_data)
        ng = data['data']['parameter']['Total number of energy groups']  
        nmat = data['data']['parameter']['Total number of materials'] 
        npc = data['data']['parameter']['Total number of pin cells']
        na = data['data']['parameter']['Total number of assemblies']  
        apc = data['data']['parameter']['Total number of active pin cells'] 
        napc = data['data']['parameter']['Total number of active pin cells']
        core = data['data']['parameter']['Core']
        if "Boundary conditions" in data['data']['parameter']:
            BC   = data['data']['parameter']['Boundary conditions']
        else:
            pass
        nx = len(core)
        for j in range(na):
            assembly.append(data['data']['Assemblies'][j]['assembly'])
        nxx = len(assembly[0])
        for j in range(npc):
            width.append(data['data']['PinCells'][j]['width'])
            regmat.append(data['data']['PinCells'][j]['mat_fill'])
            nfmesh.append(data['data']['PinCells'][j]['fine_mesh'])
        ngauss = data['data']['parameter']['Number of angular discretizations']
        Max_it = data['data']['parameter']['Maximum number of iterations']
        order =  data['data']['parameter']['The l-order Legendre polynomial']
        order = order + 1
        eps = data['data']['parameter']['Criterion of Keff convergence']
        for i in range(nmat):
            sigt.append(data['data']['materials'][i]['XSTotal'])
            nusigf.append(data['data']['materials'][i]['XSNuFission'])
            if "XSFission" in data['data']['materials'][i]:
                sigf.append(data['data']['materials'][i]['XSFission'])
            else:
                sigf.append(data['data']['materials'][i]['XSNuFission'])
            sigs.append(data['data']['materials'][i]['XSScatter Matrix'])
            chi.append(data['data']['materials'][i]['XSChi'])
        npx = len(width[0])
        totnfm = sum(nfmesh[0])*nx*nxx
        fmmid,delta = app.SlabMOC.fmm_id(assembly,core,nfmesh,width,regmat,totnfm,[npx,npc,nx,nxx,na])
        dim = ng*totnfm 
        mu,wt = app.SlabMOC.gauleg(-1.,1.,ngauss)
        p = app.SlabMOC.leg_poly(order,mu,[ngauss]) 
        d = app.SlabMOC.matrix_d(sigs,fmmid,dim,[ng,nmat,order,totnfm])
        l = app.SlabMOC.matrix_l(sigs,fmmid,dim,[ng,nmat,order,totnfm])
        u = app.SlabMOC.matrix_u(sigs,fmmid,dim,[ng,nmat,order,totnfm])
        f = app.SlabMOC.matrix_f(nusigf,chi,fmmid,dim,[ng,nmat,totnfm])  
        app.SlabMOC.title1()
        app.SlabMOC.timestamp()
        it,inter,k_eff,phi = app.SlabMOC.outer_iteration(Max_it,eps,wt,mu,d,f,u,l,p,BC,scheme2,fmmid,sigt,
                             sigf,delta,[ng,dim,totnfm,ngauss,order,nmat])
        app.SlabMOC.plot_flux(napc,delta,assembly,nfmesh,phi,fmmid,core,sigf,[dim,totnfm,nmat,ng,nx,nxx,npx,npc,na])
        sfpc = app.SlabMOC.scalarfluxpinc(nmat,ng,nfmesh,delta,assembly,phi,core,[dim,totnfm,nx,nxx,npx,npc,na])
        pf = app.SlabMOC.powerpf(napc,delta,assembly,nfmesh,phi,core,fmmid,sigf,[dim,totnfm,nmat,ng,nx,nxx,npx,npc,na])
        interval = datetime.now()-start 
        app.SlabMOC.output(str(start),str(interval),k_eff,sigt,nusigf,sigs,chi,mu,wt,width,phi,
                         eps,totnfm,it,inter,na,nx,nxx,sfpc,pf,[dim,ng,nmat,order,npx,ngauss,npc])
        print ('Total time to solution    ...........    ', interval)  
        app.SlabMOC.title2()
        del sigt,sigs,nusigf,chi,mu,wt,fmmid,nfmesh  
#=========================================================================================
#  Loading data for SN Method 1D
#========================================================================================= 
elif Methods == 'SN1D':
    with open(PATH) as json_data:
        data = json.load(json_data)
        ng = data['data']['parameter']['Total number of energy groups']  
        nmat = data['data']['parameter']['Total number of materials'] 
        npc = data['data']['parameter']['Total number of pin cells']
        na = data['data']['parameter']['Total number of assemblies']   
        napc = data['data']['parameter']['Total number of active pin cells']
        core = data['data']['parameter']['Core']
        if "Boundary conditions" in data['data']['parameter']:
            BC   = data['data']['parameter']['Boundary conditions']
        else:
            pass
        nx = len(core)
        for j in range(na):
            assembly.append(data['data']['Assemblies'][j]['assembly'])
        nxx = len(assembly[0])
        for j in range(npc):
            width.append(data['data']['PinCells'][j]['width'])
            regmat.append(data['data']['PinCells'][j]['mat_fill'])
            nfmesh.append(data['data']['PinCells'][j]['fine_mesh'])
        ngauss = data['data']['parameter']['Number of angular discretizations']
        Max_it = data['data']['parameter']['Maximum number of iterations']
        order =  data['data']['parameter']['The l-order Legendre polynomial']
        order = order + 1
        eps = data['data']['parameter']['Criterion of Keff convergence']
        for i in range(nmat):
            sigt.append(data['data']['materials'][i]['XSTotal'])
            nusigf.append(data['data']['materials'][i]['XSNuFission'])
            if "XSFission" in data['data']['materials'][i]:
                sigf.append(data['data']['materials'][i]['XSFission'])
            else:
                sigf.append(data['data']['materials'][i]['XSNuFission'])
            sigs.append(data['data']['materials'][i]['XSScatter Matrix'])
            chi.append(data['data']['materials'][i]['XSChi'])
        npx = len(width[0])
        totnfm = sum(nfmesh[0])*nx*nxx
        fmmid,delta = app.SlabSN.fmm_id(assembly,core,nfmesh,width,regmat,totnfm,[npx,npc,nx,nxx,na])
        dim = ng*totnfm 
        mu,wt = app.SlabSN.gauleg(-1.,1.,ngauss)
        p = app.SlabSN.leg_poly(order,mu,[ngauss]) 
        d = app.SlabSN.matrix_d(sigs,fmmid,dim,[ng,nmat,order,totnfm])
        l = app.SlabSN.matrix_l(sigs,fmmid,dim,[ng,nmat,order,totnfm])
        u = app.SlabSN.matrix_u(sigs,fmmid,dim,[ng,nmat,order,totnfm])
        f = app.SlabSN.matrix_f(nusigf,chi,fmmid,dim,[ng,nmat,totnfm])
        a,b = app.SlabSN.matrix_ab(dim,mu,fmmid,sigt,delta,[ng,nmat,totnfm,ngauss])
        app.SlabSN.title1()
        app.SlabSN.timestamp()
        it,inter,k_eff,phi = app.SlabSN.outer_iteration(Max_it,scheme1,eps,wt,mu,d,f,u,l,a,b,p,BC,sigt,delta,
                             nusigf,sigf,chi,fmmid,[ng,dim,totnfm,ngauss,order,nmat])
        sfpc,sf = app.SlabSN.scalarfluxpinc(nfmesh,delta,assembly,phi,sigf,fmmid,core,[dim,totnfm,nmat,ng,nx,nxx,npx,npc,na])
        app.SlabSN.plot_flux(nmat,nx,nxx,napc,delta,phi,sfpc,sf,[dim,totnfm,ng])
        pf = app.SlabSN.powerpf(totnfm,nmat,nx,nxx,napc,sfpc,sf,[ng])
        interval = datetime.now()-start
        app.SlabSN.output(str(start),str(interval),k_eff,sigt,nusigf,sigs,chi,mu,wt,width,phi,
                         eps,totnfm,it,inter,na,nx,nxx,sfpc,pf,[dim,ng,nmat,order,npx,ngauss,npc])
        print ('Total time to solution    ...........    ', interval) 
        app.SlabSN.title2()
        del sigt,sigs,nusigf,chi,mu,wt,fmmid,nfmesh
#=========================================================================================
#  Loading data for CP Method 1D
#========================================================================================= 
elif Methods == 'CP1D':
    if Geometry_type in ['Slab Geometry']:
        with open(PATH) as json_data:
            data = json.load(json_data)
            ng = data['data']['parameter']['Total number of energy groups']  
            nmat = data['data']['parameter']['Total number of materials'] 
            npc = data['data']['parameter']['Total number of pin cells']
            na = data['data']['parameter']['Total number of assemblies']   
            napc = data['data']['parameter']['Total number of active pin cells']
            core = data['data']['parameter']['Core']
            if "Boundary conditions" in data['data']['parameter']:
                BC   = data['data']['parameter']['Boundary conditions']
            else:
                pass
            nx = len(core)
            for j in range(na):
                assembly.append(data['data']['Assemblies'][j]['assembly'])
            nxx = len(assembly[0])
            for j in range(npc):
                width.append(data['data']['PinCells'][j]['width'])
                regmat.append(data['data']['PinCells'][j]['mat_fill'])
                nfmesh.append(data['data']['PinCells'][j]['fine_mesh'])
            ngauss = data['data']['parameter']['Number of angular discretizations']
            Max_it = data['data']['parameter']['Maximum number of iterations']
            order =  data['data']['parameter']['The l-order Legendre polynomial']
            order = order + 1
            eps = data['data']['parameter']['Criterion of Keff convergence']
            for i in range(nmat):
                sigtt.append(data['data']['materials'][i]['XSTotal'])
                nusigf.append(data['data']['materials'][i]['XSNuFission'])
                if "XSFission" in data['data']['materials'][i]:
                    sigf.append(data['data']['materials'][i]['XSFission'])
                else:
                    sigf.append(data['data']['materials'][i]['XSNuFission'])
                sigss.append(data['data']['materials'][i]['XSScatter Matrix'])
                chi.append(data['data']['materials'][i]['XSChi'])
            npx = len(width[0])
            totnfm = sum(nfmesh[0])*nx*nxx
            app.CP1D.title1(1)
            fmmid = app.CP1D.fmm_id(assembly,core,nfmesh,width,regmat,totnfm,[npx,npc,nx,nxx,na])
            delta,ray = app.CP1D.volume_cyl_sph_sla(1,totnfm,nfmesh,width,core,assembly,[npx,npc,na,nx,nxx])
            dim = ng*totnfm 
            if BC in ['Vacuum']:
                albedo = [0,0]
            elif BC in ['Reflective']:
                albedo = [1,1]
            elif BC in ['Vacuum Reflective']:
                albedo = [0,1]
            elif BC in ['Reflective Vacuum']:
                albedo = [1,0] 
            matrix_i = app.CP1D.matrix_matrix_i(dim) 
            sigs,sigt = app.CP1D.transport_corr(sigss,sigtt,[ng,nmat,order])
            pij = app.CP1D.pij_f(delta,albedo,sigt,fmmid,dim,[ng,totnfm,nmat])
            phi_guess = app.CP1D.flux_guess(nusigf,delta,fmmid,[dim,ng,totnfm,nmat])
            d = app.CP1D.matrix_d(sigs,fmmid,dim,[ng,totnfm,nmat])
            c = app.CP1D.matrix_c(matrix_i,d,pij,[dim])
            ainv = app.CP1D.matinv(c,[dim])
            w = app.CP1D.matrix_w(ainv,pij,[dim])
            l = app.CP1D.matrix_l(sigs,fmmid,dim,[ng,totnfm,nmat])
            u = app.CP1D.matrix_u(sigs,fmmid,dim,[ng,totnfm,nmat])
            f = app.CP1D.matrix_f(nusigf,chi,fmmid,dim,[ng,totnfm,nmat])
            a,b = app.CP1D.matrix_ab(matrix_i,l,w,u,f,[dim])
            app.CP1D.timestamp()
            iter,eval,phi = app.CP1D.aleig(a,b,eps,phi_guess,ng,totnfm,Max_it,[dim])
            app.CP1D.plot_sla_cyl_sph(1,napc,delta,assembly,nfmesh,phi,fmmid,core,sigf,[dim,totnfm,nmat,ng,nx,nxx,npx,npc,na])
            app.CP1D.powerpf(napc,delta,assembly,nfmesh,phi,core,fmmid,sigf,[dim,totnfm,nmat,ng,nx,nxx,npx,npc,na])
            sfpc = app.CP1D.scalarfluxpinc(nmat,ng,nfmesh,delta,assembly,phi,core,[dim,totnfm,nx,nxx,npx,npc,na])
            pf = app.CP1D.powerpf(napc,delta,assembly,nfmesh,phi,core,fmmid,sigf,[dim,totnfm,nmat,ng,nx,nxx,npx,npc,na])
            interval = datetime.now()-start
            app.CP1D.output(str(start),str(interval),1/eval,sigt,nusigf,sigs,chi,width,phi,eps,
                                totnfm,order,iter,iter,1,na,nx,nxx,sfpc,pf,[ng,nmat,npx,npc,dim])
            print ('Total time to solution    ...........    ', interval)  
            app.CP1D.title2()
    elif Geometry_type in ['Cylindrical Geometry']: 
        with open(PATH) as json_data:
            data = json.load(json_data)
            ng = data['data']['parameter']['Total number of energy groups']  
            nmat = data['data']['parameter']['Total number of materials'] 
            npc = data['data']['parameter']['Total number of pin cells']
            na = data['data']['parameter']['Total number of assemblies']   
            napc = data['data']['parameter']['Total number of active pin cells']
            core = data['data']['parameter']['Core']
            if "Boundary conditions" in data['data']['parameter']:
                BC   = data['data']['parameter']['Boundary conditions']
            else:
                pass
            nx = len(core)
            for j in range(na):
                assembly.append(data['data']['Assemblies'][j]['assembly'])
            nxx = len(assembly[0])
            for j in range(npc):
                width.append(data['data']['PinCells'][j]['ray'])
                regmat.append(data['data']['PinCells'][j]['mat_fill'])
                nfmesh.append(data['data']['PinCells'][j]['fine_mesh'])
            ngauss = data['data']['parameter']['Number of angular discretizations']
            Max_it = data['data']['parameter']['Maximum number of iterations']
            order =  data['data']['parameter']['The l-order Legendre polynomial']
            order = order + 1
            eps = data['data']['parameter']['Criterion of Keff convergence']
            for i in range(nmat):
                sigtt.append(data['data']['materials'][i]['XSTotal'])
                nusigf.append(data['data']['materials'][i]['XSNuFission'])
                if "XSFission" in data['data']['materials'][i]:
                    sigf.append(data['data']['materials'][i]['XSFission'])
                else:
                    sigf.append(data['data']['materials'][i]['XSNuFission'])
                sigss.append(data['data']['materials'][i]['XSScatter Matrix'])
                chi.append(data['data']['materials'][i]['XSChi'])
            npx = len(width[0])
            totnfm = sum(nfmesh[0])*nx*nxx
        app.CP1D.title1(2)
        app.CP1D.timestamp() 
        dim = ng*totnfm
        if BC in ['Vacuum']:
            albedo = 0
        elif BC in ['Reflective']:
            albedo = 1
        else:
            albedo = 2
        fmmid = app.CP1D.fmm_id(assembly,core,nfmesh,width,regmat,totnfm,[npx,npc,nx,nxx,na])
        vol,ray = app.CP1D.volume_cyl_sph_sla(2,totnfm,nfmesh,width,core,assembly,[npx,npc,na,nx,nxx])
        matrix_i = app.CP1D.matrix_matrix_i(dim) 
        track1,track3 = app.CP1D.tracking_f(0,ray,[totnfm])
        sigs,sigt = app.CP1D.transport_corr(sigss,sigtt,[ng,nmat,order])
        pij_table = app.CP1D.pij_cyl_sph(0,sigt,fmmid,track1,track3,vol,ray,albedo,dim,[ng,totnfm,nmat])
        phi_guess = app.CP1D.flux_guess(nusigf,vol,fmmid,dim,[ng,totnfm,nmat])
        d = app.CP1D.matrix_d(sigs,fmmid,dim,[ng,totnfm,nmat])
        c = app.CP1D.matrix_c(matrix_i,d,pij_table,[dim])
        ainv = app.CP1D.matinv(c,[dim])
        w = app.CP1D.matrix_w(ainv,pij_table,[dim])
        l = app.CP1D.matrix_l(sigs,fmmid,dim,[ng,totnfm,nmat])
        u = app.CP1D.matrix_u(sigs,fmmid,dim,[ng,totnfm,nmat])
        f = app.CP1D.matrix_f(nusigf,chi,fmmid,dim,[ng,totnfm,nmat])
        a,b = app.CP1D.matrix_ab(matrix_i,l,w,u,f,[dim]) 
        iter,eval,phi = app.CP1D.aleig(a,b,eps,phi_guess,ng,totnfm,Max_it,[dim])   
        app.CP1D.title2() 
        interval = datetime.now()-start 
        app.CP1D.plot_sla_cyl_sph(0,napc,ray,assembly,nfmesh,phi,fmmid,core,sigf,[dim,totnfm,nmat,ng,nx,nxx,npx,npc,na])
        sfpc = app.CP1D.scalarfluxpinc(nmat,ng,nfmesh,ray,assembly,phi,core,[dim,totnfm,nx,nxx,npx,npc,na])
        pf = app.CP1D.powerpf(napc,ray,assembly,nfmesh,phi,core,fmmid,sigf,[dim,totnfm,nmat,ng,nx,nxx,npx,npc,na])
        app.CP1D.output(str(start),str(interval),1/eval,sigt,nusigf,sigs,chi,width,phi,eps,
                                totnfm,order,iter,iter,2,na,nx,nxx,sfpc,pf,[ng,nmat,npx,npc,dim])
        print ('Total time to solution    ...........    ', interval)
    elif Geometry_type in ['Spherical Geometry']: 
        with open(PATH) as json_data:
            data = json.load(json_data)
            ng = data['data']['parameter']['Total number of energy groups']  
            nmat = data['data']['parameter']['Total number of materials'] 
            npc = data['data']['parameter']['Total number of pin cells']
            na = data['data']['parameter']['Total number of assemblies']   
            napc = data['data']['parameter']['Total number of active pin cells']
            core = data['data']['parameter']['Core']
            if "Boundary conditions" in data['data']['parameter']:
                BC   = data['data']['parameter']['Boundary conditions']
            else:
                pass
            nx = len(core)
            for j in range(na):
                assembly.append(data['data']['Assemblies'][j]['assembly'])
            nxx = len(assembly[0])
            for j in range(npc):
                width.append(data['data']['PinCells'][j]['ray'])
                regmat.append(data['data']['PinCells'][j]['mat_fill'])
                nfmesh.append(data['data']['PinCells'][j]['fine_mesh'])
            ngauss = data['data']['parameter']['Number of angular discretizations']
            Max_it = data['data']['parameter']['Maximum number of iterations']
            order =  data['data']['parameter']['The l-order Legendre polynomial']
            order = order + 1
            eps = data['data']['parameter']['Criterion of Keff convergence']
            for i in range(nmat):
                sigtt.append(data['data']['materials'][i]['XSTotal'])
                nusigf.append(data['data']['materials'][i]['XSNuFission'])
                if "XSFission" in data['data']['materials'][i]:
                    sigf.append(data['data']['materials'][i]['XSFission'])
                else:
                    sigf.append(data['data']['materials'][i]['XSNuFission'])
                sigss.append(data['data']['materials'][i]['XSScatter Matrix'])
                chi.append(data['data']['materials'][i]['XSChi'])
            npx = len(width[0])
            totnfm = sum(nfmesh[0])*nx*nxx
        app.CP1D.title1(3)
        app.CP1D.timestamp() 
        dim = ng*totnfm
        if BC in ['Vacuum']:
            albedo = 0
        elif BC in ['Reflective']:
            albedo = 1
        else:
            albedo = 2
        fmmid = app.CP1D.fmm_id(assembly,core,nfmesh,width,regmat,totnfm,[npx,npc,nx,nxx,na])
        vol,ray = app.CP1D.volume_cyl_sph_sla(3,totnfm,nfmesh,width,core,assembly,[npx,npc,na,nx,nxx])
        matrix_i = app.CP1D.matrix_matrix_i(dim) 
        track1,track3 = app.CP1D.tracking_f(1,ray,[totnfm])
        sigs,sigt = app.CP1D.transport_corr(sigss,sigtt,[ng,nmat,order])
        pij_table = app.CP1D.pij_cyl_sph(1,sigt,fmmid,track1,track3,vol,ray,albedo,dim,[ng,totnfm,nmat])
        phi_guess = app.CP1D.flux_guess(nusigf,vol,fmmid,dim,[ng,totnfm,nmat])
        d = app.CP1D.matrix_d(sigs,fmmid,dim,[ng,totnfm,nmat])
        c = app.CP1D.matrix_c(matrix_i,d,pij_table,[dim])
        ainv = app.CP1D.matinv(c,[dim])
        w = app.CP1D.matrix_w(ainv,pij_table,[dim])
        l = app.CP1D.matrix_l(sigs,fmmid,dim,[ng,totnfm,nmat])
        u = app.CP1D.matrix_u(sigs,fmmid,dim,[ng,totnfm,nmat])
        f = app.CP1D.matrix_f(nusigf,chi,fmmid,dim,[ng,totnfm,nmat])
        a,b = app.CP1D.matrix_ab(matrix_i,l,w,u,f,[dim]) 
        iter,eval,phi = app.CP1D.aleig(a,b,eps,phi_guess,ng,totnfm,Max_it,[dim])   
        app.CP1D.title2() 
        interval = datetime.now()-start 
        app.CP1D.plot_sla_cyl_sph(0,napc,ray,assembly,nfmesh,phi,fmmid,core,sigf,[dim,totnfm,nmat,ng,nx,nxx,npx,npc,na])
        sfpc = app.CP1D.scalarfluxpinc(nmat,ng,nfmesh,ray,assembly,phi,core,[dim,totnfm,nx,nxx,npx,npc,na])
        pf = app.CP1D.powerpf(napc,ray,assembly,nfmesh,phi,core,fmmid,sigf,[dim,totnfm,nmat,ng,nx,nxx,npx,npc,na])
        app.CP1D.output(str(start),str(interval),1/eval,sigt,nusigf,sigs,chi,width,phi,eps,
                                totnfm,order,iter,iter,3,na,nx,nxx,sfpc,pf,[ng,nmat,npx,npc,dim])
        print ('Total time to solution    ...........    ', interval)  

elif Methods == 'SN2D':
#=========================================================================================
#  Loading data for SN Method 2D
#========================================================================================= 
    regmat     = []
    nfmesh_xy  = []
    width_x    = []
    width_y    = []
    assembly   = []
    core       = []
    with open(PATH) as json_data:
        data = json.load(json_data)
        ng = data['data']['parameter']['Total number of energy groups']
        nmat = data['data']['parameter']['Total number of materials']
        nx = len(np.amax(data['data']['parameter']['Core'],axis=0))
        ny = len(np.amax(data['data']['parameter']['Core'],axis=1))
        npc = data['data']['parameter']['Total number of pin cells']  # Number of Pin Cells
        na = data['data']['parameter']['Total number of assemblies']  # Number of assemblies
        core = data['data']['parameter']['Core']
        napc = data['data']['parameter']['Total number of active pin cells']
        if "Boundary conditions" in data['data']['parameter']:
            BC2D   = data['data']['parameter']['Boundary conditions']
            if BC2D[0] == 'v':
                bc1 = 'Vacuum Right'
            else:
                bc1 = 'Reflective Right'
            if BC2D[1] == 'v':
                bc2 = 'Vacuum Left'
            else:
                bc2 = 'Reflective Left'
            if BC2D[2] == 'v':
                bc3 = 'Vacuum Bottom'
            else:
                bc3 = 'Reflective Bottom'
            if BC2D[3] == 'v':
                bc4 = 'Vacuum Top'
            else:
                bc4 = 'Reflective Top'
        else:
            pass
        for j in range(na):
            assembly.append(data['data']['Assemblies'][j]['assembly'])
        nxx = len(assembly[0])
        nyy = len(assembly[0])

        for j in range(npc):
            width_x.append(data['data']['PinCells'][j]['width_x'])
            width_y.append(data['data']['PinCells'][j]['width_y'])
            regmat.append(data['data']['PinCells'][j]['mat_fill'])
            nfmesh_xy.append(data['data']['PinCells'][j]['fine_mesh'])

        nord = data['data']['parameter']['Number of angular discretizations']
        Max_it = data['data']['parameter']['Maximum number of iterations']
        order =  data['data']['parameter']['The l-order Legendre polynomial']
        order = order + 1
        eps = data['data']['parameter']['Criterion of Keff convergence']
        for j in range(nmat):
            sigt.append(data['data']['materials'][j]['XSTotal'])
            nusigf.append(data['data']['materials'][j]['XSNuFission'])
            if "XSFission" in data['data']['materials'][j]:
                sigf.append(data['data']['materials'][j]['XSFission'])
            else:
                sigf.append(data['data']['materials'][j]['XSNuFission'])
            sigs.append(data['data']['materials'][j]['XSScatter Matrix'])
            chi.append(data['data']['materials'][j]['XSChi'])
        app.SN2D.title1()
        
        npx = len(width_x[0])
        npy = len(width_y[0])
        tnfm_x = sum(np.amax(nfmesh_xy[0],axis=0))*nx*nxx
        tnfm_y = sum(np.amax(nfmesh_xy[0],axis=1))*ny*nyy
        fmmid2d = app.SN2D.fmm_id2d(core,assembly,nfmesh_xy,regmat,tnfm_x,tnfm_y,[npx,npy,npc,nx,ny,nyy,nxx,na])
        xfm,yfm = app.SN2D.mesh_2d(tnfm_x,tnfm_y,core,assembly,nfmesh_xy,width_x,width_y,[nx,ny,npc,npx,npy,nxx,nyy,na])
        mu,eta,psi,w = app.SN2D.quadrature_set(nord)
        app.SN2D.timestamp()
        phi_ij,k_eff,inter,exter = app.SN2D.eigenvalues(nord,eps,Max_it,fmmid2d,mu,eta,psi,w,xfm,yfm,sigt,nusigf,sigf,
                                   bc1,bc2,bc3,bc4,chi,sigs,nfmesh_xy,[tnfm_x,tnfm_y,nmat,order,ng,ny,nx])
        app.SN2D.plot(napc,xfm,yfm,core,assembly,nusigf,sigf,chi,fmmid2d,nfmesh_xy,phi_ij,
                      [tnfm_x,tnfm_y,ny,nx,npx,npy,npc,ng,nmat,na,nxx,nyy])
        interval = datetime.now()-start 
        sfpc = app.SN2D.scalarfluxpinc(nmat,xfm,yfm,core,assembly,nfmesh_xy,phi_ij,[tnfm_y,tnfm_x,ng,nx,ny,npx,npy,npc,na,nxx,nyy])
        pf,pd = app.SN2D.pinpower(napc,xfm,yfm,core,assembly,nfmesh_xy,phi_ij,fmmid2d,nusigf,[tnfm_y,tnfm_x,nmat,ng,nx,ny,npx,npy,npc,na,nxx,nyy])
        a = app.SN2D.assemblypower(nx,ny,nxx,nyy,pf)
        app.SN2D.output(str(start),str(interval),k_eff,sigt,nusigf,sigs,chi,mu,eta,psi,w,
                         width_x,width_y,phi_ij,eps,nord,inter,exter,na,sfpc,pf,[tnfm_x,tnfm_y,ng,nmat,order,npx,npy,npc,nx*nxx,ny*nyy])      
        print ('Total time to solution   ..............  ', interval)
        app.SN2D.title2()




