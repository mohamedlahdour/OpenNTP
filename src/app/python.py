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
start = datetime.now()
sigtt  = []
sigss  = []
sigt   = []
nusigf = []
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
Comp_Method   = "Matrix Inversion"
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
#==========================================================================================
#  Loading data for MOC 1D
# ========================================================================================= 
if Methods == 'MOC':
    with open(filename) as json_data:
        data = json.load(json_data)
        ng = data['data']['parameter']['Total number of energy groups']  
        nmat = data['data']['parameter']['Total number of Materials']
        nregion = data['data']['parameter']['Total number of regions'] 
        regmat = data['data']['parameter']['Which material goes in each region']
        dcell = data['data']['parameter']['Size for each material per [cm]']
        nfmesh = data['data']['parameter']['Number of fine meshes']
        ngauss = data['data']['parameter']['Number of Angular Discretization']
        order = data['data']['parameter']['The l-order Legendre polonomial']
        order = order + 1
        Max_it = data['data']['parameter']['Maximum Iteration']
        eps = data['data']['parameter']['Epsilon Keff']
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
  
        print '  Total time to solution                     ......................................                     ', interval  
        SlabMOC.title2()  
        SlabMOC.output(str(start),BC,str(interval),k_eff,sigt,nusigf,sigs,chi,mu,wt,
                       dcell,phi,eps,totnfm,it,inter,[dim,ng,nmat,order,nregion,ngauss])
        del sigt,sigs,nusigf,chi,mu,wt,fmmid
        SlabMOC.plot_flux(delta,phi,[totnfm,dim])  
#=========================================================================================
#  Loading data for SN Method 1D
#========================================================================================= 
elif Methods == 'Sɴ':
    with open(filename) as json_data:
        data = json.load(json_data)
        ng = data['data']['parameter']['Total number of energy groups']  
        nmat = data['data']['parameter']['Total number of Materials']
        nregion = data['data']['parameter']['Total number of regions'] 
        regmat = data['data']['parameter']['Which material goes in each region']
        dcell = data['data']['parameter']['Size for each material per [cm]']
        nfmesh = data['data']['parameter']['Number of fine meshes']
        ngauss = data['data']['parameter']['Number of Angular Discretization']
        order = data['data']['parameter']['The l-order Legendre polonomial']
        order = order + 1
        Max_it = data['data']['parameter']['Maximum Iteration']
        eps = data['data']['parameter']['Epsilon Keff']
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
        phi_ni,phi_li = SlabSN.flux_guess(dim,totnfm,nusigf,dcell,wt,p,[ng,nmat,nregion,ngauss,order])
        mp,mm = SlabSN.matrix_mp(ng,nmat,totnfm,ngauss,a,b,[dim])
        a = SlabSN.matrix_a(ng,nmat,totnfm,mp,mm,[dim,ngauss])
        SlabSN.title1()
        SlabSN.timestamp()
        if Comp_Method in ["Gauss Seidel, Liebmann with Successive Over-Relaxation"]:
            iter,inter,k_eff,phi = SlabSN.gsl_sor(Max_it,eps,wt,mu,d,f,u,l,p,a,sigt,phi_ni,
                                   phi_li,delta,[ng,dim,totnfm,ngauss,order,nmat])
        if Comp_Method in ["Matrix Inversion"]:  
            n = totnfm*ngauss
            ainv = SlabSN.inverse(a,[n])
            it,inter,k_eff,phi = SlabSN.eigenvalues(Max_it,eps,wt,mu,d,f,u,l,p,ainv,sigt,phi_ni,phi_li,delta,
                                                   [ng,dim,totnfm,ngauss,order,nmat])
        if Comp_Method in ["Conjugate gradient"]:
            iter,inter,k_eff,phi = SlabSN.conjgrad(Max_it,eps,wt,mu,d,f,u,l,p,a,sigt,phi_ni,phi_li,delta,
                                                    [ng,dim,totnfm,ngauss,order,nmat])
            #del fmmid,nfmesh
            #it,inter,k_eff,phi = SlabSN.outer_iteration(Max_it,scheme1,eps,wt,mu,d,f,u,l,a,b,p,BC,sigt,
            #                                   flux_ni,flux_li,delta,[ng,dim,totnfm,ngauss,order,nmat])
        interval = datetime.now()-start        
        print '  Total time to solution      ........................     ', interval   
        SlabSN.title2()
        #SlabSN.output(str(start),BC,str(interval),k_eff,sigt,nusigf,sigs,chi,mu,wt,
        #              dcell,phi,eps,totnfm,it,inter,[dim,ng,nmat,order,nregion,ngauss])
        #del sigt,sigs,nusigf,chi,mu,wt
        SlabSN.plot_flux(delta,phi,[totnfm,dim])
#=========================================================================================
#  Loading data for CP Method 1D
#========================================================================================= 
elif Methods == 'CP':
    if Geometry_type in ['Slab Geometry']:
        with open(filename) as json_data:
            data = json.load(json_data)
            ng = data['data']['parameter']['Total number of energy groups']        
            nmat = data['data']['parameter']['Total number of Materials']
            nregion = data['data']['parameter']['Total number of regions'] 
            regmat = data['data']['parameter']['Which material goes in each region']
            dcell = data['data']['parameter']['Size for each material per [cm]']
            nfmesh = data['data']['parameter']['Number of fine meshes']
            Max_it = data['data']['parameter']['Maximum Iteration']
            order =  data['data']['parameter']['The l-order Legendre polonomial']
            order = order + 1
            eps = data['data']['parameter']['Epsilon Keff']
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
        print '  Total time to solution      ........................     ', interval    
        CP1D.title2() 
        CP1D.output(str(start),albedo,str(interval),1./eval,sigt,nusigf,sigs,chi,dcell,
                                          phi,eps,totnfm,order,iter,[dim,ng,nmat,nregion])
        del sigt,nusigf,sigs,chi,dcell,fmmid,l,d,c,ainv,w,u,f
        CP1D.plot_sla_cyl_sph(vol,1,phi,[totnfm,dim])

    if Geometry_type in ['Cylindrical Geometry']: 
        with open(filename) as json_data:
            data = json.load(json_data)
            ng = data['data']['parameter']['Total number of energy groups'] 
            nmat = data['data']['parameter']['Total number of Materials']
            nregion = data['data']['parameter']['Total number of regions'] 
            regmat = data['data']['parameter']['Which material goes in each region']
            rayon = data['data']['parameter']['Ray for each region per [cm]']
            nfmesh = data['data']['parameter']['Number of fine meshes']
            Max_it = data['data']['parameter']['Maximum Iteration']
            order =  data['data']['parameter']['The l-order Legendre polonomial']
            order = order + 1
            eps = data['data']['parameter']['Epsilon Keff']
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
        print vol
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
            #print a
        CP1D.title1(2)
        CP1D.timestamp() 
        iter,eval,phi = CP1D.aleig(a,b,eps,phi_guess,ng,totnfm,Max_it,[dim])   
        CP1D.title2() 
        CP1D.plot_sla_cyl_sph(ray,0,phi,[totnfm,dim])

    if Geometry_type in ['Spherical Geometry']: 
        with open(filename) as json_data:
            data = json.load(json_data)
            ng = data['data']['parameter']['Total number of energy groups'] 
            nmat = data['data']['parameter']['Total number of Materials']
            nregion = data['data']['parameter']['Total number of regions'] 
            regmat = data['data']['parameter']['Which material goes in each region']
            rayon = data['data']['parameter']['Ray for each region per [cm]']
            nfmesh = data['data']['parameter']['Number of fine meshes']
            Max_it = data['data']['parameter']['Maximum Iteration']
            order =  data['data']['parameter']['The l-order Legendre polonomial']
            order = order + 1
            eps = data['data']['parameter']['Epsilon Keff']
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
        CP1D.title2() 
        CP1D.plot_sla_cyl_sph(ray,0,phi,[totnfm,dim])
elif Methods == 'SN':
#=========================================================================================
#  Loading data for SN Method 2D
#========================================================================================= 
    regmat     = []
    nfmesh_xy  = []
    with open(filename) as json_data:
        data = json.load(json_data)
        ng = data['data']['parameter']['Total number of energy groups'] 
        nmat = data['data']['parameter']['Total number of Materials']
        nx = data['data']['parameter']['Total number of X regions']
        ny = data['data']['parameter']['Total number of Y regions'] 
        xcm = data['data']['parameter']['X region thickness per [cm]']
        ycm = data['data']['parameter']['Y region thickness per [cm]']
        for j in range(ny):
            regmat.append(data['data']['parameter']['Which material goes in each cell'][j])
            nfmesh_xy.append(data['data']['parameter']['XY number of fine meshes in each cell'][j])
        nord = data['data']['parameter']['Number of Angular Discretization']
        Max_it = data['data']['parameter']['Maximum Iteration']
        order =  data['data']['parameter']['The l-order Legendre polonomial']
        order = order + 1
        eps = data['data']['parameter']['Epsilon Keff']
        for j in range(nmat):
            sigtt.append(data['data']['materials'][j]['XSTotal'])
            nusigf.append(data['data']['materials'][j]['XSNuFission'])
            sigss.append(data['data']['materials'][j]['XSScatter Matrix'])
            chi.append(data['data']['materials'][j]['XSChi']) 
        SlabSN2D.title1() 
        tnfm_x = sum(np.amax(nfmesh_xy,axis=0)) # x
        tnfm_y = sum(np.amax(nfmesh_xy,axis=1)) # y
        fmmid2d = SlabSN2D.fmm_id(nfmesh_xy,regmat,tnfm_y,tnfm_x,[nx,ny])
        xfm,yfm = SlabSN2D.mesh2d(tnfm_x,tnfm_y,nfmesh_xy,xcm,ycm,[nx,ny])
        mu,eta,psi,w = SlabSN2D.quadrature_set(nord)
        #phi_lmij = SlabSN2D.guessflux(xfm,yfm,nusigf,chi,fmmid2d,order,[tnfm_x,tnfm_y,nmat,ng])
        #qs_ij = SlabSN2D.scattering_source(nord,fmmid2d,mu,eta,w,sigss,phi_lmij,[tnfm_x,tnfm_y,nmat,order,ng])
        #dsnew,qf_ij = SlabSN2D.fission_source(nord,fmmid2d,xfm,yfm,mu,w,nusigf,chi,phi_lmij,1.0,[tnfm_x,tnfm_y,nmat,order,ng])
        #print (qf_ij)
        #qs_ij = SlabSN2D.scattering_source(ngauss,fmmid2d,mu,w,sigss,phi_lmij,[tnfm_x,tnfm_y,nmat,order,ng]) 
        #dsnew,qf_ij = SlabSN2D.fission_source(ngauss,fmmid2d,xfm,yfm,mu,w,nusigf,
        #                           chi,phi_lmij,1.0,[tnfm_x,tnfm_y,nmat,order,ng])
        #plm = SlabSN2D.plm(1,1,mu)
        #print (plm)
        #rlm = SlabSN2D.rlm(1,1,mu,eta)
        #print (rlm, eta,mu)
        SlabSN2D.timestamp()
        phi_ij,k_eff,inter,exter = SlabSN2D.eigenvalues(nord,fmmid2d,mu,eta,w,xfm,yfm,sigtt,nusigf,
                 bc1,bc2,bc3,bc4,chi,sigss,nfmesh_xy,[tnfm_x,tnfm_y,nmat,order,ng,ny,nx])

        SlabSN2D.plot_flux(xfm,yfm,fmmid2d,nfmesh_xy,phi_ij,[tnfm_x,tnfm_y,ny,nx,ng])
        #p = SlabSN2D.powerdensity(nfmesh_xy,regmat,fmmid2d,nusigf,chi,xfm,
        #                          yfm,phi_ij,[nx,ny,tnfm_y,tnfm_x,nmat,ng])
        interval = datetime.now()-start  
        SlabSN2D.output(start,interval,k_eff,sigtt,nusigf,sigss,chi,mu,eta,psi,w,
                        xcm,ycm,phi_ij,eps,nord,inter,exter,[tnfm_y,tnfm_x,ng,nmat,order,nx,ny])      
        print 'Total time to solution   ..............  ', interval 
        SlabSN2D.title2() 
        


        #a = SlabSN2D.plot_plm(order-1,order-1)

