#! /usr/bin/env python
#! -*- coding:utf-8 -*-cd
import CP1D
import SlabSN
import SlabMOC
import SlabSN2D
import numpy as np
import json
from datetime import datetime
import subprocess
start  = datetime.now()
sigtt  = []
sigss  = []
sigt   = []
nusigf = []
sigf   = []
sigs   = []
chi    = []
vol    = []
filename      = open('app/link/script.py', "r" ).read()
Methods       = open('app/link/script00.py', "r" ).read()
Geometry_type = open('app/link/script01.py', "r" ).read()
BC            = open('app/link/script02.py', "r" ).read()
scheme1       = open('app/link/script03.py', "r" ).read()
scheme2       = open('app/link/script04.py', "r" ).read()
bc1           = open('app/link/script06.py', "r" ).read()
bc2           = open('app/link/script07.py', "r" ).read()
bc3           = open('app/link/script08.py', "r" ).read()
bc4           = open('app/link/script09.py', "r" ).read()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
#==========================================================================================
#  Loading data for MOC 1D
# ========================================================================================= 
if Methods == 'MOC1D':
    with open(filename) as json_data:
        data = json.load(json_data)
        ng = data['data']['parameter']['Total number of energy groups']  
        nmat = data['data']['parameter']['Total number of Materials']
        nregion = data['data']['parameter']['Total number of regions'] 
        regmat = data['data']['parameter']['Which material fills each region']
        dcell = data['data']['parameter']['Size of each region [cm]']
        nfmesh = data['data']['parameter']['Number of fine meshes']
        ngauss = data['data']['parameter']['Number of Angular Discretization']
        order = data['data']['parameter']['The l-order Legendre polynomial']
        order = order + 1
        Max_it = data['data']['parameter']['Maximum Number of Iterations']
        eps = data['data']['parameter']['Criterion of Keff convergence']
        for i in range(nmat):
            sigt.append(data['data']['materials'][i]['XSTotal'])
            nusigf.append(data['data']['materials'][i]['XSNuFission'])
            sigs.append(data['data']['materials'][i]['XSScatter Matrix'])
            chi.append(data['data']['materials'][i]['XSChi'])
        totnfm = sum(nfmesh)
        dim = ng*totnfm 
        mu,wt = SlabMOC.gauleg(-1.,1.,ngauss)
        p = SlabMOC.leg_poly(order,mu,[ngauss]) 
        fmmid = SlabMOC.fmm_id(totnfm,nfmesh,regmat,[nregion])
        delta = SlabMOC.delta_f(totnfm,nfmesh,dcell,[nregion])
        d = SlabMOC.matrix_d(sigs,fmmid,dim,[ng,nmat,order,totnfm])
        l = SlabMOC.matrix_l(sigs,fmmid,dim,[ng,nmat,order,totnfm])
        u = SlabMOC.matrix_u(sigs,fmmid,dim,[ng,nmat,order,totnfm])
        f = SlabMOC.matrix_f(nusigf,chi,fmmid,dim,[ng,nmat,totnfm])
        flux_ni,flux_li = SlabMOC.flux_guess(dim,nregion,nusigf,dcell,p,
                                            [ng,nmat,ngauss,order,totnfm]) 
        del nfmesh            
        SlabMOC.title1()
        SlabMOC.timestamp()
        it,inter,k_eff,phi = SlabMOC.outer_iteration(Max_it,eps,wt,mu,d,f,u,l,p,BC,scheme2,fmmid,sigt,
                             flux_ni,flux_li,delta,[ng,dim,totnfm,ngauss,order,nmat])
        interval = datetime.now()-start  
        print 'Total time to solution   ..............  ', interval  
        SlabMOC.title2()
        SlabMOC.output(str(start),BC,str(interval),k_eff,sigt,nusigf,sigs,chi,mu,wt,
                       dcell,phi,eps,totnfm,it,inter,[dim,ng,nmat,order,nregion,ngauss])
        del sigt,sigs,nusigf,chi,mu,wt,fmmid
        SlabMOC.plot_flux(delta,phi,[totnfm,dim])   
#=========================================================================================
#  Loading data for SN Method 1D
#========================================================================================= 
elif Methods == 'SN1D':
    with open(filename) as json_data:
        data = json.load(json_data)
        ng = data['data']['parameter']['Total number of energy groups']  
        nmat = data['data']['parameter']['Total number of Materials']
        nregion = data['data']['parameter']['Total number of regions'] 
        regmat = data['data']['parameter']['Which material fills each region']
        dcell = data['data']['parameter']['Size of each region [cm]']
        nfmesh = data['data']['parameter']['Number of fine meshes']
        ngauss = data['data']['parameter']['Number of Angular Discretization']
        order = data['data']['parameter']['The l-order Legendre polynomial']
        order = order + 1
        Max_it = data['data']['parameter']['Maximum Number of Iterations']
        eps = data['data']['parameter']['Criterion of Keff convergence']
        for i in range(nmat):
            sigt.append(data['data']['materials'][i]['XSTotal'])
            nusigf.append(data['data']['materials'][i]['XSNuFission'])
            sigs.append(data['data']['materials'][i]['XSScatter Matrix'])
            chi.append(data['data']['materials'][i]['XSChi'])
        totnfm = sum(nfmesh)
        dim = ng*totnfm 
        mu,wt = SlabSN.gauleg(-1.,1.,ngauss)
        p = SlabSN.leg_poly(order,mu,[ngauss]) 
        fmmid = SlabSN.fmm_id(totnfm,nfmesh,regmat,[nregion])
        delta = SlabSN.delta_f(totnfm,nfmesh,dcell,[nregion])
        d = SlabSN.matrix_d(sigs,fmmid,dim,[ng,nmat,order,totnfm])
        l = SlabSN.matrix_l(sigs,fmmid,dim,[ng,nmat,order,totnfm])
        u = SlabSN.matrix_u(sigs,fmmid,dim,[ng,nmat,order,totnfm])
        f = SlabSN.matrix_f(nusigf,chi,fmmid,dim,[ng,nmat,totnfm])
        a,b = SlabSN.matrix_ab(dim,mu,fmmid,sigt,delta,[ng,nmat,totnfm,ngauss])
        flux_ni,flux_li = SlabSN.flux_guess(dim,nregion,nusigf,dcell,wt,p,
                                           [ng,nmat,ngauss,order,totnfm])
        del fmmid,nfmesh
        SlabSN.title1()
        SlabSN.timestamp()
        it,inter,k_eff,phi = SlabSN.outer_iteration(Max_it,scheme1,eps,wt,mu,d,f,u,l,a,b,p,BC,sigt,
                             flux_ni,flux_li,delta,[ng,dim,totnfm,ngauss,order,nmat])

        interval = datetime.now()-start  
        print 'Total time to solution   ..............  ', interval   
        SlabSN.title2()
        SlabSN.output(str(start),BC,str(interval),k_eff,sigt,nusigf,sigs,chi,mu,wt,
                      dcell,phi,eps,totnfm,it,inter,[dim,ng,nmat,order,nregion,ngauss])
        del sigt,sigs,nusigf,chi,mu,wt
        SlabSN.plot_flux(delta,phi,[totnfm,dim])
#=========================================================================================
#  Loading data for CP Method 1D
#========================================================================================= 
elif Methods == 'CP1D':
    if Geometry_type in ['Slab Geometry']:
        with open(filename) as json_data:
            data = json.load(json_data)
            ng = data['data']['parameter']['Total number of energy groups']        
            nmat = data['data']['parameter']['Total number of Materials']
            nregion = data['data']['parameter']['Total number of regions'] 
            regmat = data['data']['parameter']['Which material fills each region']
            dcell = data['data']['parameter']['Size of each region [cm]']
            nfmesh = data['data']['parameter']['Number of fine meshes']
            Max_it = data['data']['parameter']['Maximum Number of Iterations']
            order =  data['data']['parameter']['The l-order Legendre polynomial']
            order = order + 1
            eps = data['data']['parameter']['Criterion of Keff convergence']
            totnfm = sum(nfmesh)
            for i in range(nmat):
                sigtt.append(data['data']['materials'][i]['XSTotal'])
                nusigf.append(data['data']['materials'][i]['XSNuFission'])
                sigss.append(data['data']['materials'][i]['XSScatter Matrix'])
                chi.append(data['data']['materials'][i]['XSChi']) 
        if BC in ['Vacuum']:
            albedo = [0,0]
        elif BC in ['Reflective']:
            albedo = [1,1]
        elif BC in ['Vacuum Reflective']:
            albedo = [0,1]
        elif BC in ['Reflective Vacuum']:
            albedo = [1,0]      
        for i in range(nregion):
            for j in range(nfmesh[i]):
                vol.append(dcell[i]/nfmesh[i])
        dim = ng*totnfm
        fmmid = CP1D.fmm_id(nfmesh,regmat,totnfm,[nregion])
        matrix_i = CP1D.matrix_matrix_i(dim) 
        sigs,sigt = CP1D.transport_corr(sigss,sigtt,[ng,nmat,order])    
        pij = CP1D.pij_f(vol,albedo,sigt,fmmid,dim,[ng,totnfm,nmat])
        phi_guess = CP1D.flux_guess(nusigf,vol,fmmid,[dim,ng,totnfm,nmat])
        d = CP1D.matrix_d(sigs,fmmid,dim,[ng,totnfm,nmat])
        c = CP1D.matrix_c(matrix_i,d,pij,[dim])
        ainv = CP1D.matinv(c,[dim])
        w = CP1D.matrix_w(ainv,pij,[dim])
        l = CP1D.matrix_l(sigs,fmmid,dim,[ng,totnfm,nmat])
        u = CP1D.matrix_u(sigs,fmmid,dim,[ng,totnfm,nmat])
        f = CP1D.matrix_f(nusigf,chi,fmmid,dim,[ng,totnfm,nmat])
        a,b = CP1D.matrix_ab(matrix_i,l,w,u,f,[dim])
        CP1D.title1(1) 
        CP1D.timestamp()
        iter,eval,phi = CP1D.aleig(a,b,eps,phi_guess,ng,totnfm,Max_it,[dim])  
        interval = datetime.now()-start  
        print 'Total time to solution   ..............  ', interval    
        CP1D.title2() 
        CP1D.output(str(start),albedo,str(interval),1./eval,sigt,nusigf,sigs,chi,dcell,
                                          phi,eps,totnfm,order,iter,1,[dim,ng,nmat,nregion])
        del sigt,nusigf,sigs,chi,dcell,fmmid,l,d,c,ainv,w,u,f
        CP1D.plot_sla_cyl_sph(vol,1,phi,[totnfm,dim])

    if Geometry_type in ['Cylindrical Geometry']: 
        with open(filename) as json_data:
            data = json.load(json_data)
            ng = data['data']['parameter']['Total number of energy groups'] 
            nmat = data['data']['parameter']['Total number of Materials']
            nregion = data['data']['parameter']['Total number of regions'] 
            regmat = data['data']['parameter']['Which material fills each region']
            rayon = data['data']['parameter']['Ray for each region per [cm]']
            nfmesh = data['data']['parameter']['Number of fine meshes']
            Max_it = data['data']['parameter']['Maximum Number of Iterations']
            order =  data['data']['parameter']['The l-order Legendre polynomial']
            order = order + 1
            eps = data['data']['parameter']['Criterion of Keff convergence']
            totnfm = sum(nfmesh)
            for i in range(nmat):
                sigtt.append(data['data']['materials'][i]['XSTotal'])
                nusigf.append(data['data']['materials'][i]['XSNuFission'])
                sigss.append(data['data']['materials'][i]['XSScatter Matrix'])
                chi.append(data['data']['materials'][i]['XSChi'])  
            
        dim = ng*totnfm
        if BC in ['Vacuum']:
            albedo = 0
        elif BC in ['Reflective']:
            albedo = 1
        elif BC in ['Vacuum Reflective']:
            albedo = 2
        elif BC in ['Reflective Vacuum']:
            albedo = 2
        fmmid = CP1D.fmm_id(nfmesh,regmat,totnfm ,[nregion])
        ray,vol = CP1D.vol_ray(0,nfmesh,rayon,totnfm,[nregion])
        matrix_i = CP1D.matrix_matrix_i(dim) 
        track1,track3 = CP1D.tracking_f(0,ray,[totnfm])
        sigs,sigt = CP1D.transport_corr(sigss,sigtt,[ng,nmat,order])
        pij_table = CP1D.pij_cyl_sph(0,sigt,fmmid,track1,track3,vol,ray,albedo,dim,[ng,totnfm,nmat])
        #print pij_table
        phi_guess = CP1D.flux_guess(nusigf,vol,fmmid,dim,[ng,totnfm,nmat])
        #print phi_guess
        d = CP1D.matrix_d(sigs,fmmid,dim,[ng,totnfm,nmat])
        c = CP1D.matrix_c(matrix_i,d,pij_table,[dim])
        ainv = CP1D.matinv(c,[dim])
        w = CP1D.matrix_w(ainv,pij_table,[dim])
        l = CP1D.matrix_l(sigs,fmmid,dim,[ng,totnfm,nmat])
        u = CP1D.matrix_u(sigs,fmmid,dim,[ng,totnfm,nmat])
        f = CP1D.matrix_f(nusigf,chi,fmmid,dim,[ng,totnfm,nmat])
        a,b = CP1D.matrix_ab(matrix_i,l,w,u,f,[dim])
        interval = datetime.now()-start  
        print 'Total time to solution   ..............  ', interval  
        CP1D.title1(2)
        CP1D.timestamp() 
        iter,eval,phi = CP1D.aleig(a,b,eps,phi_guess,ng,totnfm,Max_it,[dim])   
        CP1D.title2() 
        CP1D.output(str(start),[0,0],str(interval),1./eval,sigt,nusigf,sigs,chi,rayon,
                                          phi,eps,totnfm,order,iter,2,[dim,ng,nmat,nregion])
        CP1D.plot_sla_cyl_sph(ray,0,phi,[totnfm,dim])

    if Geometry_type in ['Spherical Geometry']: 
        with open(filename) as json_data:
            data = json.load(json_data)
            ng = data['data']['parameter']['Total number of energy groups'] 
            nmat = data['data']['parameter']['Total number of Materials']
            nregion = data['data']['parameter']['Total number of regions'] 
            regmat = data['data']['parameter']['Which material fills each region']
            rayon = data['data']['parameter']['Ray for each region per [cm]']
            nfmesh = data['data']['parameter']['Number of fine meshes']
            Max_it = data['data']['parameter']['Maximum Number of Iterations']
            order =  data['data']['parameter']['The l-order Legendre polynomial']
            order = order + 1
            eps = data['data']['parameter']['Criterion of Keff convergence']
            totnfm = sum(nfmesh)
            for i in range(nmat):
                sigtt.append(data['data']['materials'][i]['XSTotal'])
                nusigf.append(data['data']['materials'][i]['XSNuFission'])
                sigss.append(data['data']['materials'][i]['XSScatter Matrix'])
                chi.append(data['data']['materials'][i]['XSChi'])  
        dim = ng*totnfm
        if BC in ['Vacuum']:
            albedo = 0
        elif BC in ['Reflective']:
            albedo = 1
        elif BC in ['Vacuum Reflective']:
            albedo = 2
        elif BC in ['Reflective Vacuum']:
            albedo = 2
        ray,vol = CP1D.vol_ray(1,nfmesh,rayon,totnfm,[nregion])
        fmmid = CP1D.fmm_id(nfmesh,regmat,totnfm,[nregion])
        matrix_i = CP1D.matrix_matrix_i(dim)
        track1,track3 = CP1D.tracking_f(1,ray,[totnfm])
        sigs,sigt = CP1D.transport_corr(sigss,sigtt,[ng,nmat,order])
        pij_table = CP1D.pij_cyl_sph(1,sigt,fmmid,track1,track3,vol,ray,albedo,dim,[ng,totnfm,nmat])
        phi_guess = CP1D.flux_guess(nusigf,vol,fmmid,[dim,ng,totnfm,nmat])
        d = CP1D.matrix_d(sigs,fmmid,dim,[ng,totnfm,nmat])
        c = CP1D.matrix_c(matrix_i,d,pij_table,[dim])
        ainv = CP1D.matinv(c,[dim])
        w = CP1D.matrix_w(ainv,pij_table,[dim])
        l = CP1D.matrix_l(sigs,fmmid,dim,[ng,totnfm,nmat])
        u = CP1D.matrix_u(sigs,fmmid,dim,[ng,totnfm,nmat])
        f = CP1D.matrix_f(nusigf,chi,fmmid,dim,[ng,totnfm,nmat])
        a,b = CP1D.matrix_ab(matrix_i,l,w,u,f,[dim])
        CP1D.title1(3)
        CP1D.timestamp()
        iter,eval,phi = CP1D.aleig(a,b,eps,phi_guess,ng,totnfm,Max_it,[dim])
        interval = datetime.now()-start
        print 'Total time to solution   ..............  ', interval
        CP1D.title2()
        CP1D.output(str(start),[0,0],str(interval),1./eval,sigt,nusigf,sigs,chi,rayon,
                                          phi,eps,totnfm,order,iter,3,[dim,ng,nmat,nregion])
        CP1D.plot_sla_cyl_sph(ray,0,phi,[totnfm,dim])
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
    with open(filename) as json_data:
        data = json.load(json_data)
        ng = data['data']['parameter']['Total number of energy groups']
        nmat = data['data']['parameter']['Total number of materials']
        nx = len(np.amax(data['data']['parameter']['Core'],axis=0))
        ny = len(np.amax(data['data']['parameter']['Core'],axis=1))
        npc = data['data']['parameter']['Total number of pin cells']  # Number of Pin Cells
        na = data['data']['parameter']['Total number of assemblies']  # Number of assemblies
        core = data['data']['parameter']['Core']
        napc = data['data']['parameter']['Total number of active pin cells']

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
                sigf=[[0]*ng]*nmat
            sigs.append(data['data']['materials'][j]['XSScatter Matrix'])
            chi.append(data['data']['materials'][j]['XSChi'])
        SlabSN2D.title1()
        
        npx = len(width_x[0])
        npy = len(width_y[0])
        tnfm_x = sum(np.amax(nfmesh_xy[0],axis=0))*nx*nxx
        tnfm_y = sum(np.amax(nfmesh_xy[0],axis=1))*ny*nyy
        assembl = SlabSN2D.convert(core,assembly,[nx,ny,na,nxx,nyy])
        fmmid2d = SlabSN2D.fmm_id(assembl,nfmesh_xy,regmat,tnfm_x,tnfm_y,[npx,npy,npc,nx*nxx,ny*nyy])
        xfm,yfm = SlabSN2D.mesh2d(tnfm_x,tnfm_y,assembl,nfmesh_xy,width_x,width_y,[nx*nxx,ny*nyy,npc,npx,npy])
        mu,eta,psi,w = SlabSN2D.quadrature_set(nord)

        SlabSN2D.timestamp()
        phi_ij,k_eff,inter,exter = SlabSN2D.eigenvalues(nord,eps,Max_it,fmmid2d,mu,eta,psi,w,xfm,yfm,sigt,nusigf,sigf,
                                    bc1,bc2,bc3,bc4,chi,sigs,nfmesh_xy,[tnfm_x,tnfm_y,nmat,order,ng,ny*nyy,nx*nxx])
        SlabSN2D.plot(napc,xfm,yfm,assembl,nusigf,sigf,chi,fmmid2d,nfmesh_xy,phi_ij,
                      [tnfm_x,tnfm_y,ny*nyy,nx*nxx,npx,npy,npc,ng,nmat])
        interval = datetime.now()-start 
        sfpc = SlabSN2D.scalarfluxpinc(nmat,xfm,yfm,assembl,nfmesh_xy,phi_ij,[tnfm_x,tnfm_y,ng,nx*nxx,ny*nyy,npx,npy,npc])
        pf,pd = SlabSN2D.powerpf(xfm,yfm,napc,assembl,nfmesh_xy,phi_ij,fmmid2d,nusigf,[tnfm_x,tnfm_y,nmat,ng,nx*nxx,ny*nyy,npx,npy,npc])
        SlabSN2D.output(str(start),str(interval),k_eff,sigt,nusigf,sigs,chi,mu,eta,psi,w,
                         width_x,width_y,phi_ij,eps,nord,inter,exter,na,sfpc,pf,[tnfm_x,tnfm_y,ng,nmat,order,npx,npy,npc,nx*nxx,ny*nyy])      
        print 'Total time to solution   ..............  ', interval
        SlabSN2D.title2()
