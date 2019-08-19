!BL1
    subroutine quadrature_set(N,mup, etap, psip,pwi)
       !Calculate the quadrature sets in an octant and return
       !complete angle and weight set over all
       implicit none
       integer(kind=8), intent(in) :: N
       real(kind=8), dimension(N/2) :: mmu
       real(kind=8), dimension(N*(N+2)/8), intent(out) :: mup, etap, psip,pwi
       real(kind=8), parameter :: pi = 3.141592653589793
       real(kind=8) :: mu1,C
       integer(kind=8) :: i,j,k,ip,M
       character(3) :: CR
       character(100) :: ligne
       open (99,file='app/Output/Ordinates.h')
       open (98,file='app/data/weights.h')
       M = N*(N+2)/8
       select case (N)
          case (2)
             mu1=0.33333333; CR = '1  '
          case (4)
             mu1=0.12251480; CR = '3  ' 
          case (6)
             mu1=0.07109447; CR = '6  '
          case (8)
             mu1=0.04761903; CR = '10 '
          case (10)
             mu1=0.03584310; CR = '15 '
          case (12)
             mu1=0.02796615; CR = '21 '
          case (14)
             mu1=0.02310250; CR = '28 '
          case (16)
             mu1=0.01931398; CR = '36 '
          case (18)
             mu1=0.01692067; CR = '45 '
          case (20)
             mu1=0.01455253; CR = '55 '
          case default 
             print * , " SN METHOD 2D: ORDER NOT AVAILABLE FOR N."
             stop
       end select
       
           mmu(1) = sqrt(mu1)
           C = 2.0D0*(1.0D0-3.0D0*mu1)/DBLE(N-2)
       if ( N > 2) then
           do i = 2,N/2
              mmu(i) = sqrt(mu1 + C*DBLE(i-1))
           enddo
       endif
      ip=0
      k=1
      do j = N/2,1,-1
         do i=1,k
         ip=ip+1 
         mup(ip)  = mmu(j)
         enddo
         k=k+1
      end do

     if (N==2) then
         etap(1)  = mmu(1)
         psip(1) = mmu(1)
     else
      ip=0
         do i=1,2
         ip=ip+1 
         etap(ip)  = mmu(1)
         psip(ip) = mmu(i)
         enddo
         ip=ip+1
         etap(ip)  = mmu(2)
         psip(ip) = mmu(1)
      if (N>4) then
         k=3
      10 do j=1,k
         ip=ip+1
         etap(ip)  = mmu(j)
         psip(ip)  =  mmu(k-j+1)
         enddo
         k=k+1
         if (k>N/2) then
             goto 11
         else
             go to 10
         endif
      endif
   11 endif    
       do
          ligne = ""
          read(98,"(a)")ligne 
          k = index( ligne , CR)
          ! si k n'est pas 0 c'est que ligne contient "CR" 
          if ( k /= 0 ) exit
       enddo
       do i=1,N*(N+2)/8 
          read(98,*) pwi(i)
       enddo

       do i = 1,N*(N+2)/8  
          write(99,3000) mup(i),etap(i),psip(i),pwi(i)
       enddo
       3000 format(4(f12.7)) 
    end subroutine quadrature_set
!BL2
    subroutine fmm_id(fmmid2D,nfmesh_xy,RegMat,nx,ny,TNFM_y,TNFM_x)
       implicit none
       integer(kind=4), intent(in) :: nx,ny,TNFM_y,TNFM_x
       integer(kind=4), dimension(ny,nx), intent(in) :: RegMat,nfmesh_xy
       integer(kind=4), dimension(TNFM_x,TNFM_y), intent(out) :: fmmid2D
       integer(kind=4), dimension(TNFM_x,TNFM_y) :: fmmid
       integer(kind=4) :: i,j,m,n,nn,mm
       !-- fine mesh material id x and y directions  
       mm=1
       do i = 1,nx
          nn = 1
          do j = 1,ny
             n = nfmesh_xy(j,1); m = nfmesh_xy(1,i)
                 fmmid(mm:m+mm-1,nn:n+nn-1) = RegMat(j,i)
                 nn=nn+n
          enddo
          mm=mm+m
       enddo
       fmmid2D = fmmid(1:TNFM_x,TNFM_y:1:-1)
       !write(*,'(16i2)') fmmid2D
    end subroutine fmm_id
!BL3    
    subroutine Mesh2D(nx,ny,TNFM_x,TNFM_y,nfmesh_xy,xcm,ycm,xfm,yfm)
       implicit none
       integer(kind=4), intent(in) :: nx,ny,TNFM_x,TNFM_y
       integer(kind=4), dimension(ny,nx),    intent(in) :: nfmesh_xy
       real(kind=8),    dimension(nx),       intent(in)  :: xcm
       real(kind=8),    dimension(ny),       intent(in)  :: ycm
       real(kind=8),    dimension(TNFM_x), intent(out) :: xfm
       real(kind=8),    dimension(TNFM_y), intent(out) :: yfm
       integer(kind=4) :: i,j,k,n
       !>  nx            Total number of Coarse mesh along x axis.
       !>  ny            Total number of Coarse mesh along y axis.
       !>  TNFM_x      Total number of fine mesh along x axis.
       !>  TNFM_y      Total number of fine mesh along y axis.
       !>  xfm           Number of fine mesh per x coarse mesh.
       !>  yfm           Number of fine mesh per y coarse mesh.
       !>  xcm           Coarse mesh edges along x axis.
       !>  ycm           Coarse mesh edges along y axis.
       !>  Size each mesh in the x direction
       n = 1
       do i = 1,nx
          do k = 1,nfmesh_xy(1,i)
             xfm(n) = xcm(i)/nfmesh_xy(1,i)
             n = n+1
          enddo
       enddo
       !-- Size each mesh in the y direction
       n = 1
       do j = 1,ny
          do k = 1,nfmesh_xy(j,1)
             yfm(n) = ycm(j)/nfmesh_xy(j,1)
             n = n+1
          enddo
       enddo
    end subroutine Mesh2D
!BL4
    subroutine GuessFlux(phi_lmij,TNFM_x,TNFM_y,xfm,yfm,NusigF,Chi,fmmid2D,&
                         order,Nmat,ng)
        implicit none
        integer(kind=4), intent(in) :: TNFM_x,TNFM_y,order,Nmat,ng
        integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
        real(kind=8), dimension(Nmat,ng), intent(in) :: NusigF,Chi
        real(kind=8), dimension(TNFM_x), intent(in) :: xfm
        real(kind=8), dimension(TNFM_y), intent(in) :: yfm
        real(kind=8), dimension(TNFM_x,TNFM_y,order,order+1,ng), intent(out) :: phi_lmij
        real(kind=8) :: val
        integer(kind=4) :: i,j,k
        phi_lmij = 0.0D0
        do k = 1,ng      
           do j = 1,TNFM_y
              do i = 1,TNFM_x
                 val = xfm(i)*yfm(j)*sum(Chi(:,k)*sum(NusigF(:,k)))
                 phi_lmij(i,j,1,1,k) = 1.0d0/val
              enddo
           enddo
        enddo
    end subroutine GuessFlux
!BL5
    subroutine Scattering_Source(TNFM_x,TNFM_y,Nmat,order,ng,NORD,&
                                 fmmid2D,mu,eta,w,SigS,phi_lmij,Qs_ij)
        implicit none
        integer(kind=4), intent(in) :: TNFM_x,TNFM_y,Nmat,order,ng,NORD
        real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
        real(kind=8), dimension(NORD*(NORD+2)/8), intent(in) :: mu,eta,w
        integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
        real(kind=8), dimension(TNFM_x,TNFM_y,order,order+1,ng), intent(in) :: phi_lmij
        real(kind=8), dimension(TNFM_x,TNFM_y,ng), intent(out) :: Qs_ij
        ! Variables locales
        real(kind=8), dimension(TNFM_x,TNFM_y,order,order+1,ng) :: Qs_lmij
        real(kind=8), dimension(TNFM_x,TNFM_y,NORD*(NORD+2)/2,ng) :: Qs_nij
        real(kind=8), dimension(TNFM_x,TNFM_y,ng) :: Qs
        real(kind=8), parameter :: pi = 3.141592653589793
        real(kind=8) :: Rlm,s1,s2
        integer(kind=4) :: i,j,k,l,m,n,n1,i1,numord,mm
        numord=NORD*(NORD+2)/8
       
        !>--------------- 
        Qs_nij = 0.0d0; Qs = 0.0d0; Qs_lmij=0

        do k=1,ng
        do j = 1,TNFM_y
        do i = 1,TNFM_x
           do l = 0,order-1
              mm=1
              do m = -l,l
                 if (mod(l+m,2) .eq. 0) then
                 Qs_lmij(i,j,l+1,mm,k) = sum(SigS(fmmid2D(i,j),l+1,:,k)&
                                                 *phi_lmij(i,j,l+1,mm,:))
                 mm=mm+1
                 endif
              enddo
           enddo
        enddo
        enddo
        enddo 

       
        !>--------------- 
        do k = 1,ng
        do j = 1,TNFM_y
        do i = 1,TNFM_x
        n1 = 1
        do i1 = 1,4
           do n = 1,numord
                  s2 = 0.0d0
                  mm=1
                  do l = 0,order-1   
                     s1 = 0.0d0
                     do m=-l,l
                        if (mod(l+m,2) .eq. 0) then
                        s1 = s1 + Qs_lmij(i,j,l+1,mm,k)*Rlm(l,m,mu(n),eta(n))
                        mm=mm+1
                        endif
                     enddo
                     s2 = s2 + (2.0d0*l+1.0d0)*s1
                  enddo
                  Qs_nij(i,j,n1,k) = Qs_nij(i,j,n1,k) + s2
           n1 = n1+1             
           enddo  
        enddo
        enddo
        enddo
        enddo
    
        !>---------------
        do k = 1,ng
        do j = 1,TNFM_y
           do i = 1,TNFM_x
             n1=1
             do i1 = 1,4
               do n  = 1,numord
                  Qs(i,j,k) = Qs(i,j,k) + 0.25D0*w(n)*Qs_nij(i,j,n1,k) !> end
                  n1=n1+1
               enddo
             enddo
           enddo
        enddo
        enddo 
        Qs_ij = Qs     
    end subroutine Scattering_Source
!BL6
    subroutine Fission_Source(TNFM_x,TNFM_y,Nmat,order,ng,NORD,fmmid2D,&
                              xfm,yfm,dsnew,mu,w,NusigF,Chi,phi_lmij,k_eff,Qf_ij)
        implicit none
        integer(kind=4), intent(in) :: TNFM_x,TNFM_y,Nmat,order,ng,NORD
        real(kind=8), intent(in) :: k_eff 
        real(kind=8), dimension(Nmat,ng), intent(in) :: NusigF,Chi
        real(kind=8), dimension(NORD*(NORD+2)/8), intent(in) :: mu,w
        integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
        real(kind=8), dimension(TNFM_x,TNFM_y,order,order+1,ng), intent(in) :: phi_lmij
        real(kind=8), dimension(TNFM_x), intent(in) :: xfm
        real(kind=8), dimension(TNFM_y), intent(in) :: yfm
        real(kind=8), dimension(TNFM_x,TNFM_y,ng), intent(out) :: Qf_ij
        real(kind=8), intent(out) :: dsnew
        ! Varaibles locales
        real(kind=8), dimension(TNFM_x,TNFM_y,NORD*(NORD+2)/2,ng) :: Qf_nij
        real(kind=8), dimension(TNFM_x,TNFM_y,order,order+1,ng) :: phi_11ij, Qf_lmij
        real(kind=8), dimension(TNFM_x*TNFM_y*ng) :: V
        real(kind=8), dimension(TNFM_x,TNFM_y,ng) :: Qf
        real(kind=8), parameter :: PI = 3.141592653589793                                   
        real(kind=8) :: som,plm
        integer(kind=4) :: i,j,k,l,m,n,n1,i1,numord
        numord = NORD*(NORD+2)/8
        phi_11ij = 0.0d0
        phi_11ij(:,:,1,1,:) = phi_lmij(:,:,1,1,:)
        Qf_nij = 0.0d0; Qf = 0.0d0
        !>-------------
        do k=1,ng
        do j = 1,TNFM_y
        do i = 1,TNFM_x
                 Qf_lmij(i,j,1,1,k) = Chi(fmmid2D(i,j),k)*sum(NusigF(fmmid2D(i,j),:)&
                                          *phi_11ij(i,j,1,1,:))/k_eff 
        enddo
        enddo
        enddo
        !>--------------- 
        do k = 1,ng
        do j = 1,TNFM_y
        do i = 1,TNFM_x
           do n = 1,NORD*(NORD+2)/2
                  Qf_nij(i,j,n,k) = Qf_nij(i,j,n,k) + Qf_lmij(i,j,1,1,k)
           enddo  
        enddo
        enddo
        enddo
        !>---------------
        m=1
        do k = 1,ng
        do j = 1,TNFM_y
        do i = 1,TNFM_x
             n1=1
             do i1 = 1,4
               do n  = 1,numord
                  Qf(i,j,k) = Qf(i,j,k) + 0.25D0*w(n)*Qf_nij(i,j,n1,k) !> end
                  n1=n1+1
               enddo
             enddo
             V(m) = Qf(i,j,k)
             m = m+1
        enddo
        enddo
        enddo 
        Qf_ij = Qf
        !V=V/sum(V)
        dsnew = sum(V)
        
    end subroutine Fission_Source
!BL7
    subroutine Eigenvalues(TNFM_x,TNFM_y,Nmat,order,ng,N,ny,nx,&
                           fmmid2D,mu,eta,w,xfm,yfm,SigT,NusigF,&
                           BC1,BC2,BC3,BC4,Chi,SigS,phi_ij,nfmesh_xy,k_eff,inter,exter)
    implicit none
    ! Variables Globales
    CHARACTER(50), intent(in) :: BC1,BC2,BC3,BC4
    integer(kind=4), intent(in) :: TNFM_x,TNFM_y,Nmat,order,ng,N,ny,nx
    integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
    real(kind=8), dimension(N*(N+2)/8), intent(in) :: mu,eta,w
    real(kind=8), dimension(TNFM_x), intent(in) :: xfm
    real(kind=8), dimension(TNFM_y), intent(in) :: yfm
    real(kind=8), dimension(Nmat,ng), intent(in) :: SigT,NusigF,Chi
    real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
    integer(kind=4), dimension(ny,nx), intent(in) :: nfmesh_xy
    real(kind=8), dimension(TNFM_x,TNFM_y,ng), intent(out) :: phi_ij
    real(kind=8), intent(out) :: k_eff
    integer(kind=4), intent(out) :: inter,exter
    ! Variables Locales
    real(kind=8), dimension(TNFM_x,TNFM_y,N*(N+2)/2,ng) :: phi_ijn
    real(kind=8), dimension(TNFM_x,TNFM_y,order,order+1,ng) :: phi_lmij
    real(kind=8), dimension(TNFM_x,TNFM_y,ng) :: Qf_ij,Qs_ij,Q_ij
    real(kind=8), dimension(TNFM_x,TNFM_y,ng) :: phi_ij0,phi_ij1
    real(kind=8), dimension(TNFM_x*TNFM_y*ng) :: V,V1
    real(kind=8) :: err_k_eff, err_phi,eps,k_eff0,err_flux
    real(kind=8) :: dsnew, dsold,Rlm
    integer(kind=4) :: i,j,k,l,m,n1,i1,n2,mm
    !> -----------------
    call GuessFlux(phi_lmij,TNFM_x,TNFM_y,xfm,yfm,NusigF,Chi,fmmid2D,&
                     order,Nmat,ng)   
    call NormalizeFlux(TNFM_x,TNFM_y,Nmat,ng,NusigF,Chi,fmmid2D,xfm,yfm,phi_lmij) 
    phi_ij1   = phi_lmij(:,:,1,1,:)
    phi_ij0   = phi_lmij(:,:,1,1,:)
    inter     = 0
    err_k_eff = 1.0D0
    err_phi   = 1.0D0
    exter     = 0 
    k_eff     = 1.0D0 ! Given K_eff 
    k_eff0    = 1.0D0
    eps       = 1.0E-6
    call Fission_Source(TNFM_x,TNFM_y,Nmat,order,ng,N,fmmid2D,&
                        xfm,yfm,dsnew,mu,w,NusigF,Chi,phi_lmij,k_eff,Qf_ij)
    dsold = dsnew
    do while (err_k_eff >= eps .and. err_phi >= 1.0E-7)  ! Staring extern Iteration

       call Fission_Source(TNFM_x,TNFM_y,Nmat,order,ng,N,fmmid2D,&
                              xfm,yfm,dsnew,mu,w,NusigF,Chi,phi_lmij,k_eff,Qf_ij)
       err_flux = 1.0D0
       do while (err_flux >= 1.0E-7)  ! Starting intern Iteration
          call Scattering_Source(TNFM_x,TNFM_y,Nmat,order,ng,N,fmmid2D,&
                                 mu,eta,w,SigS,phi_lmij,Qs_ij)
          Q_ij = Qf_ij + Qs_ij
          call Sweep2D(BC1,BC2,BC3,BC4,TNFM_x,TNFM_y,N,ng,&
                            Nmat,SigT,Q_ij,mu,eta,w,xfm,yfm,fmmid2D,phi_ij,phi_ijn)  
          !phi_lmij(:,:,1,1,:) = phi_ij(:,:,:)
          phi_lmij = 0.0D0
          mm=1
          do l=0,order-1
          do m=-l,l
             if (mod(l+m,2) .eq. 0) then
          n1=1
          do i1=1,4
          do n2=1,N*(N+2)/8
             phi_lmij(:,:,l+1,mm,:)=phi_lmij(:,:,l+1,mm,:)+0.25D0*w(n2)*phi_ijn(:,:,n1,:)*Rlm(l,m,mu(n2),eta(n2))
             n1=n1+1
          enddo
          enddo
          mm=mm+1
          endif
          enddo
          enddo
       
 

          err_flux = maxval(abs(phi_ij-phi_ij0)/phi_ij0) ! Condition sur le flux angulaire
          phi_ij0 = phi_ij
          inter = inter + 1
       enddo ! ending intern Iteration
       err_phi =  maxval(abs(phi_ij-phi_ij1)/phi_ij1)  ! Condition sur le flux scalaire
       if (exter==0) go to 100    
       k_eff =  dsnew/dsold*k_eff
       err_k_eff =  abs(k_eff-k_eff0)/k_eff
       write(*,2000)exter,k_eff,err_k_eff
       100 continue
       dsold    =  dsnew
       k_eff0   = k_eff
       phi_ij1  = phi_ij
       exter    = exter + 1
    enddo ! ending extern Iteration
       phi_ij=phi_ij/maxval(phi_ij) 
    2000 format("Iteration",i4,":",1x,"===>",1x,"keff =",F9.6,1x,"===>",1x,"res =",e10.3)
    end subroutine Eigenvalues
!BL8
    subroutine plot_flux(TNFM_x,TNFM_y,ny,nx,ng,xfm,yfm,fmmid2D,nfmesh_xy,phi_ij)
       implicit none
       integer(kind=4), intent(in) :: TNFM_x,TNFM_y,ny,nx,ng
       real(kind=8), dimension(TNFM_x), intent(in) :: xfm
       real(kind=8), dimension(TNFM_y), intent(in) :: yfm
       integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
       real(kind=8), dimension(TNFM_x,TNFM_y,ng), intent(in) :: phi_ij
       integer(kind=4), dimension(ny,nx), intent(in) :: nfmesh_xy
       real(kind=8), dimension(TNFM_y) :: dy
       real(kind=8), dimension(TNFM_x) :: dx
       real(kind=8) :: r
       integer(kind=4) :: i,j,n,k,m,k1
       open (14,file='app/Output/flux_sn2D.h') 
       m = 1 ;r=0.0d0
       do j=1,ny
             do n=1,nfmesh_xy(j,1)
                dy(m) = sum(yfm(1:m))
                m=m+1
             enddo
       enddo
       m = 1 
       do i=1,nx
          do n=1,nfmesh_xy(1,i)
             dx(m) = sum(xfm(1:m))
             m=m+1
          enddo
       enddo

       write(14,'(f15.8,10000f15.8)') r,(dx(i), i=1,TNFM_x) 
       do j = 1,TNFM_y
          write(14,'(f15.8,10000f15.8)') dy(j),(phi_ij(i,j,1), i=1,TNFM_x)  
       enddo
       close(14)
    end subroutine plot_flux
!BL----------------------------
    subroutine Sweep2D(BC1,BC2,BC3,BC4,TNFM_x,TNFM_y,NORD,ng,&
                       Nmat,SigT,Q_ij,mu,eta,w,dx,dy,fmmid2D,phi_ij,phi_ijn)
    implicit none
    ! Variable Globale
    CHARACTER(50), intent(in) :: BC1,BC2,BC3,BC4
    integer(kind=4), intent(in) :: TNFM_x,TNFM_y,NORD,ng,Nmat
    real(kind=8), dimension(TNFM_x,TNFM_y,ng), intent(in) :: Q_ij 
    real(kind=8), dimension(NORD*(NORD+2)/8), intent(in) :: mu, eta, w
    real(kind=8), dimension(TNFM_x), intent(in) :: dx
    real(kind=8), dimension(TNFM_y), intent(in) :: dy
    real(kind=8), dimension(Nmat,ng),  intent(in) :: SigT
    integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
    real(kind=8), dimension(TNFM_x,TNFM_y,NORD*(NORD+2)/2,ng), intent(out) :: phi_ijn
    real(kind=8), dimension(TNFM_x,TNFM_y,ng), intent(out) :: phi_ij
    ! Variable locale
    real(kind=8), parameter :: PI = 3.141592653589793
    real(kind=8), dimension(TNFM_x,TNFM_y,NORD*(NORD+2)/2,ng) :: phi_aijn,phi_bijn
    real(kind=8), dimension(TNFM_x,TNFM_y,ng) :: phi
    real(kind=8), dimension(TNFM_x+1,TNFM_y,NORD*(NORD+2)/2) :: phix,phiax,phibx
    real(kind=8), dimension(TNFM_x,TNFM_y+1,NORD*(NORD+2)/2) :: phiy,phiay,phiby
    real(kind=8) :: ax,ay,nume,deno,beta
    integer(kind=4) :: i,j,n,k,numord,nn
    numord = NORD*(NORD+2)/2
    phi = 0.0d0
    !       -----------------------------------------------------------------        
    !        mu < 0, eta < 0,   RIGHT-TO-LEFT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    phix(TNFM_x+1,TNFM_y:1:-1,1:numord/4) = 0.0D0
    phiy(TNFM_x:1:-1,TNFM_y+1,1:numord/4) = 0.0D0 
    do k = 1,ng
    do n = 1,numord/4
       do j = TNFM_y,1,-1  
          ay = 2.0D0*abs(eta(n))/dy(j)       
          do i = TNFM_x,1,-1     
             ax = 2.0D0*abs(mu(n))/dx(i)
             deno = SigT(fmmid2D(i,j),k) + ax + ay 
             nume = Q_ij(i,j,k) + ax*phix(i+1,j,n) + ay*phiy(i,j+1,n)
             phi_ijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0D0 !> end
             phix(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i+1,j,n)
             !> simple negative flux fixup         
             if (phix(i,j,n)<0) phix(i,j,n) = 0.0D0 !> end
             !>  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k) !> end
          enddo
          do i = TNFM_x,1,-1         
             phiy(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j+1,n)
             !> simple negative flux fixup         
             if (phiy(i,j,n)<0) phiy(i,j,n) = 0.0D0 !> end
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------        
    !        mu > 0, eta < 0,   LEFT-TO-RIGHT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    phix(1,TNFM_y:1:-1,numord/4+1:numord/2) = phix(1,TNFM_y:1:-1,numord/4:1:-1)    
    phiy(1:TNFM_x,TNFM_y+1,numord/4+1:numord/2) = 0.0d0 
    do k = 1,ng
    nn=1
    do n = numord/4+1,numord/2     
       do j = TNFM_y,1,-1               
          ay = 2.0d0*abs(eta(nn))/dy(j)   
          do i = 1,TNFM_x
             ax = 2.0d0*abs(mu(nn))/dx(i)
             deno = SigT(fmmid2D(i,j),k) + ax + ay 
             nume = Q_ij(i,j,k) + ax*phix(i,j,n) + ay*phiy(i,j+1,n)
             phi_ijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0D0 !> end
             phix(i+1,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i,j,n)
             !> simple negative flux fixup         
             if (phix(i+1,j,n)<0) phix(i+1,j,n) = 0.0D0 !> end
             !> update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(nn)*phi_ijn(i,j,n,k) !> end
          enddo
          do i = 1,TNFM_x
             phiy(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j+1,n)
             !> simple negative flux fixup         
             if (phiy(i,j,n)<0) phiy(i,j,n) = 0.0d0 !> end
          enddo
       enddo
       nn=nn+1
    enddo
    enddo
    !       -----------------------------------------------------------------
    !        mu < 0, eta > 0,   RIGHT-TO-LEFT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    phix(TNFM_x+1,1:TNFM_y,3*numord/4+1:numord) = 0.0d0 
    phiy(TNFM_x:1:-1,1,3*numord/4+1:numord) = phiy(TNFM_x:1:-1,1,numord/4:1:-1)    
    do k = 1,ng
    nn=1
    do n = 3*numord/4+1,numord 
       do j = 1,TNFM_y 
          ay = 2.0d0*abs(eta(nn))/dy(j)   
          do i = TNFM_x,1,-1      
             ax = 2.0d0*abs(mu(nn))/dx(i)
             deno = SigT(fmmid2D(i,j),k) + ax + ay 
             nume = Q_ij(i,j,k) + ax*phix(i+1,j,n) + ay*phiy(i,j,n)
             phi_ijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0D0 !> end
             phix(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i+1,j,n)
             !> simple negative flux fixup         
             if (phix(i,j,n)<0) phix(i,j,n) = 0.0D0 !> end
             !>  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(nn)*phi_ijn(i,j,n,k)
          enddo
          do i = TNFM_x,1,-1         
             phiy(i,j+1,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j,n)
             !> simple negative flux fixup         
             if (phiy(i,j+1,n)<0) phiy(i,j+1,n) = 0.0D0 !> end
          enddo
       enddo
       nn=nn+1
    enddo
    enddo
    !       ----------------------------------------------------------------- 
    !        mu > 0, eta > 0,   LEFT-TO-RIGHT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    phix(1,1:TNFM_y,numord/2+1:3*numord/4) = phix(1,1:TNFM_y,numord:3*numord/4:-1) 
    phiy(1:TNFM_x,1,numord/2+1:3*numord/4) = phiy(1:TNFM_x,1,numord/2:numord/4:-1)   
    do k = 1,ng
    nn=1
    do n = numord/2+1,3*numord/4   
       do j = 1,TNFM_y      
          ay = 2.0D0*abs(eta(nn))/dy(j)
          do i = 1,TNFM_x  
             ax = 2.0D0*abs(mu(nn))/dx(i)
             deno = SigT(fmmid2D(i,j),k) + ax + ay 
             nume = Q_ij(i,j,k) + ax*phix(i,j,n) + ay*phiy(i,j,n)
             phi_ijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0D0 !> end
             phix(i+1,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i,j,n)
             !> simple negative flux fixup         
             if (phix(i+1,j,n)<0) phix(i+1,j,n) = 0.0D0 !> end
             !>  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(nn)*phi_ijn(i,j,n,k) !> end
          enddo
          do i = 1,TNFM_x
             phiy(i,j+1,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j,n)
             !> simple negative flux fixup         
             if (phiy(i,j+1,n)<0) phiy(i,j+1,n) = 0.0D0 !> end
          enddo
       enddo
       nn=nn+1
    enddo
    enddo
    phi_ij = phi
    end subroutine Sweep2D

!BL21
    subroutine timestamp()
!      ------------------------------------------------------------------------
!      TIMESTAMP prints the current YMDHMS date as a time stamp.
!      Example:
!      31 May 2001   9:45:54.872 AM
!      Licensing:
!      This code is distributed under the GNU LGPL license.
!      Modified:
!      18 May 2013
!      Author:
!      John Burkardt
!      Parameters:
!      None
!      ------------------------------------------------------------------------
       implicit none

       character ( len = 8 ) ampm
       integer ( kind = 4 ) d
       integer ( kind = 4 ) h
       integer ( kind = 4 ) m
       integer ( kind = 4 ) mm
       character ( len = 9 ), parameter, dimension(12) :: month = (/ &
       'January  ', 'February ', 'March    ', 'April    ', &
       'May      ', 'June     ', 'July     ', 'August   ', &
       'September', 'October  ', 'November ', 'December ' /)
       integer ( kind = 4 ) n
       integer ( kind = 4 ) s
       integer ( kind = 4 ) values(8)
       integer ( kind = 4 ) y
       call date_and_time ( values = values )
       y = values(1)
       m = values(2)
       d = values(3)
       h = values(5)
       n = values(6)
       s = values(7)
       mm = values(8)
       if ( h < 12 ) then
       ampm = 'AM'
       else if ( h == 12 ) then
       if ( n == 0 .and. s == 0 ) then
       ampm = 'Noon'
       else
       ampm = 'PM'
       end if
       else
       h = h - 12
       if ( h < 12 ) then
       ampm = 'PM'
       else if ( h == 12 ) then
       if ( n == 0 .and. s == 0 ) then
       ampm = 'Midnight'
       else
       ampm = 'AM'
       end if
       end if
       end if
       write ( *, '(i6,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
       d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
       return
    end subroutine timestamp
!BL11
    subroutine title1() 
       print*, '  ____                   _____   _____ _   _ '  
       print*, ' / __ \                 |  __ \ / ____| \ | |'
       print*, '| |  | |_ __   ___ _ __ | |__) | (___ |  \| |'
       print*, '| |  | | "_ \ / _ \ "_ \|  _  / \___ \| . ` |'
       print*, '| |__| | |_) |  __/ | | | | \ \ ____) | |\  |'
       print*, ' \____/| .__/ \___|_| |_|_|  \_\_____/|_| \_|'
       print*, '       | |                                   '
       print*, '       |_|                                   '
       print*, '____________________________________________________________'
       print*, '         | The OpenRSN Neutron Transport Package            '       
       print*, 'Version  | Version Number: 1.2                              '     
       print*, 'Copyright| 2015-2019 Radiation and Nuclear system Laboratory'
       print*, '         | University Abdelmalk Essaadi, FS Tetouan, Morocco'
       print*, 'Source   | FORTRAN90 version                                ' 
       print*, 'GUI      | PyQt5                                            ' 
       print*, 'Method   | The Discrete Ordinates Method Sn                 '  
       print*, 'Dimension| Two dimensions (2D)                              '
       print*, 'Geometry | Cartesian                                        ' 
       print*, '____________________________________________________________'
       print*, ''
    end subroutine title1
!BL12
    subroutine title2()
       write(*,FMT=*)'**********************************************************'
       write(*,FMT=*)'                         Finished                         '                             
       write(*,FMT=*)'**********************************************************'  
    end subroutine title2
!BL13
    function plm(l,m,mu) 
    ! calcululation of associated Legendre functions
    !  0 ≤ m ≤ ℓ
    implicit none
    integer(kind=4), intent(in) :: l,m
    real(kind=8), intent(in) :: mu
    real(kind=8) :: x=0,plm
    if (m<0) then
       print*,'bad arguments in alf (m<0)'
       stop
    elseif (m>l) then
       print*,'bad arguments in alf (m>l)'
       stop
    elseif (abs(mu)>1.) then
       print*,'bad arguments in alf (abs(mu)>1.)'
       stop
    endif
          if (l==0 .and. m==0) then
              x = 1.
          else if (l==1 .and. m==0) then
              x = mu
          else if (l==1 .and. m==1) then
              x = -sqrt(1. - mu**2)
          else if (l==2 .and. m==0) then
              x = 0.5*(3.*mu**2- 1.)
          else if (l==2 .and. m==1) then
              x = -3*mu*sqrt(1.0-mu**2)
          else if (l==2 .and. m==2) then
              x = 3*(1.0-mu**2)
          else if (l==3 .and. m==0) then
              x = 0.5*sqrt(5*mu**3-3*mu)
          else if (l==3 .and. m==1) then
              x = -1.5*(5*mu**2-1.)*sqrt(1.-mu**2)
          else if (l==3 .and. m==2) then
              x = 15*mu*(1.-mu**2)
          else if (l==3 .and. m==3) then
              x = -15*(1.-mu**2)**(3./2.)
          else if (l==4 .and. m==0) then
              x = 0.125*(35*mu**4-30*mu**2+3)
          else if (l==4 .and. m==1) then
              x = -2.5*(7*mu**3-1.)*sqrt(1.-mu**2)
          else if (l==4 .and. m==2) then
              x = 7.5*(7*mu**2-1.)*(1.-mu**2)
          else if (l==4 .and. m==3) then
              x = -105*mu*(1.-mu**2)**(3./2.)
          else if (l==4 .and. m==4) then
              x = 105*(1.-mu**2)**2
          endif 
          plm = x
    end function 
!BL14
    subroutine NormalizeFlux(TNFM_x,TNFM_y,Nmat,ng,NusigF,Chi,fmmid2D,xfm,yfm,phi) 
    implicit none
    ! Neutron scalar flux is normalized according to sum(V*Chi*NusigF*phi=1)
    integer(kind=4), intent(in) :: TNFM_x,TNFM_y,Nmat,ng
    real(kind=8), dimension(Nmat,ng), intent(in) :: NusigF,Chi
    real(kind=8), dimension(TNFM_x), intent(in) :: xfm
    real(kind=8), dimension(TNFM_y), intent(in) :: yfm
    integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
    real(kind=8), dimension(TNFM_x,TNFM_y), intent(inout) :: phi
    real(kind=8), dimension(TNFM_x*TNFM_y) :: V
    integer(kind=4) :: i,j,m
    real(kind=8) :: s
    m=1
    do j = 1,TNFM_y
       do i = 1,TNFM_x
          V(m) = xfm(i)*yfm(j)*Chi(fmmid2D(i,j),1)*NusigF(fmmid2D(i,j),1)*phi(i,j)
          m=m+1
       enddo
    enddo
    phi = phi/sum(V)

    s=0
    do j = 1,TNFM_y
       do i = 1,TNFM_x
          s = s + xfm(i)*yfm(j)*Chi(fmmid2D(i,j),1)*NusigF(fmmid2D(i,j),1)*phi(i,j)
       enddo
    enddo
    end subroutine NormalizeFlux
!BL15
    function Rlm(l,m,mu,eta) 
       implicit none
       integer(kind=4), intent(in) :: l,m
       real(kind=8), intent(in) :: mu,eta
       real(kind=8) :: Rlm,plm,phi,Im,delta,fact1,fact2
       integer(kind=4) :: i
       phi=acos(eta/sqrt(abs(1-mu**2)))
       if (m >=0) then
          Im=cos(m*phi)
       else
          Im=sin(abs(m)*phi)
       endif
       if (m==0) then
          delta = 1.0d0
       else
          delta = 0.0d0
       endif
       fact1=1
       fact2=1 
       do i=1,(l-abs(m)) 
          fact1=fact1*i
       end do
       do i=1,(l+abs(m)) 
          fact2=fact2*i
       end do
       Rlm = abs(sqrt((2.0d0-delta)*fact1/fact2)*plm(l,abs(m),mu)*Im) 
    end function
!BL16
    subroutine PowerDensity(nx,ny,TNFM_y,TNFM_x,Nmat,ng,nfmesh_xy,&
                               RegMat,fmmid2D,NusigF,Chi,xfm,yfm,phi_ij,P)
    ! CALCULATION OF POWER DENSITY DISTRIBUTION 
       implicit none
       integer(kind=4), intent(in) :: nx,ny,TNFM_y,TNFM_x,Nmat,ng
       integer(kind=4), dimension(ny,nx), intent(in) :: RegMat,nfmesh_xy
       integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
       real(kind=8), dimension(Nmat,ng), intent(in) :: NusigF,Chi
       real(kind=8), dimension(TNFM_x), intent(in) :: xfm
       real(kind=8), dimension(TNFM_y), intent(in) :: yfm
       real(kind=8), dimension(TNFM_x,TNFM_y,ng), intent(in) :: phi_ij
       real(kind=8), dimension(ny,nx), intent(out) :: P
       integer(kind=8) :: i,j,n,m,nn,mm
       mm=1
       do i = 1,nx
          nn = 1
          do j = 1,ny
             n = nfmesh_xy(j,1); m = nfmesh_xy(1,i)
             P(i,j) = sum(xfm(mm:m+mm-1)*yfm(nn:n+nn-1)*Chi(RegMat(j,i),1)*&
                      NusigF(RegMat(j,i),1)*phi_ij(mm:m+mm-1,j,1))
                 nn=nn+n
          enddo
             mm=mm+m
       enddo
       write(*,'(f8.4)') P
    end subroutine PowerDensity    
!BL17
subroutine Output(start,tm,k_eff,SigT,NusigF,SigS,Chi,mup,etap,psip,pwi,xcm,ycm,&
                 phi_ij,eps,TNFM_y,TNFM_x,ng,Nmat,order,nx,ny,N,it1,it2)
        implicit none
        integer(kind=4), intent(in) :: ng,TNFM_y,TNFM_x,Nmat,order,nx,ny,N,it1,it2
        real(kind=8), dimension(Nmat,ng), intent(in) :: SigT,NusigF,Chi
        real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
        real(kind=8), dimension(TNFM_x,TNFM_y,ng), intent(in) :: phi_ij
        real(kind=8), dimension(nx), intent(in)  :: xcm
        real(kind=8), dimension(ny), intent(in)  :: ycm
        real(kind=8), dimension(N*(N+2)/8), intent(in) :: mup,etap,psip,pwi
        CHARACTER(50), intent(in) :: start,tm
        real(kind=8), intent(in) :: eps,k_eff
        integer(kind=4) :: i,j
        open (100,file='app/Output/OUTPUT_SN2D.TXT')
        write (100, FMT=* ) '********************************************************************************'
        write (100, FMT=* ) 'ERSN, UNIVERSITY ABDELMALEK ESSAADI FACULTY OF SCIENCES - TETOUAN, MOROCCO'
        write (100, FMT=* ) 'CODE  DEVELOPED  BY  MOHAMED  LAHDOUR,  PHD  STUDENT'
        write (100, FMT=* ) 'OpenRSN:         SN  DISCRETE  ORDINATES  METHOD'
        write (100, FMT=* ) 'DIMENSION:       TWO DIMENSIONS (2D) '
        write (100, FMT=* ) 'GEOMETRY:        CARTESIAN'
        write (100, FMT=* ) 'VERSION NUMBER:  1.2'
        write (100, FMT=* ) 'VERSION DATE:    20  August  2019'
                             

        write (100,3010) 'RAN ON:          ', start,'(H:M:S)'
        write (100, FMT=* ) '********************************************************************************'
        write (100, FMT=* ) '           ----------------------------------------------------------' 
        write (100, FMT=* ) '                     INPUT  PARAMETER - VALUES  FROM  INPUT'              
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'ENERGY GROUPS NUMBER:                    ',ng
        write (100, FMT=* ) 'X REGIONS NUMBER:                        ',nx
        write (100, FMT=* ) 'Y REGIONS NUMBER:                        ',ny
        write (100, FMT=* ) 'MATERIALS NUMBER:                        ',Nmat
        write (100,3040)    'SIZE OF EACH X REGION [CM]:              ',xcm 
        write (100,3040)    'SIZE OF EACH Y REGION [CM]:              ',ycm     
        write (100, FMT=* ) 'NUMBER OF DIRECTION ALONG EACH AXIS:     ',N
        write (100, FMT=* ) 'ORDER LEGENDRE POLYNOMIAL:               ',order-1
        write (100, FMT=* ) 'TOTAL NUMBER OF X FINE MESHES:           ',TNFM_x
        write (100, FMT=* ) 'TOTAL NUMBER OF Y FINE MESHES:           ',TNFM_y
        write (100,3050)    'CONVERGENCE CRITERION of KEFF AND FLUX:  ',eps
        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '                      CALCULATION  RUN-TIME  PARAMETERS  SN' 
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'LEVEL  SYMMETRIC  GAUSSIAN  QUADRATURE  SETS: '
        write (100, FMT=* ) ''
        write (100, FMT=* ) '      N. ORDER ','         MU    ','         ETA    ','         PSI    ','     WEIGHTS '
        write (100, FMT=* ) ''
        do i = 1,N*(N+2)/8  
          write(100,3060) i,mup(i),etap(i),psip(i),pwi(i)
        enddo
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'PSEUDO  CROSS  SECTIONS  DATA: '
        write (100, FMT=* ) ''

        do i = 1,Nmat
        write (100, 3070) ' MATERIAL :', i  
        write (100, FMT=* ) ''
        write (100, FMT=* ) '        GROUP ','          TOTAL ','       ABSORPTION ',&
                            '     NU*FISSION ','     SCATTERING ','     FISSION SPECTRUM'
        write (100, FMT=* ) ''
            do j = 1,ng
            write(100,3080) j,SigT(i,j),SigT(i,j)-SigS(i,1,j,j),NusigF(i,j),SigS(i,1,j,j),Chi(i,j)
            enddo
        enddo
        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '                       NORMALIZED SCALAR  FLUX  SOLUTION' 
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'FLUXES  PER  MESH  PER  ENERGY  GROUP:'  
        write (100, FMT=* ) '' 
        write (100,3000)' M E S H ', ('     G R O U P',i,i=1,ng)
        write (100, FMT=* ) ''
        write(100,'(i14,100i13)') (i, i=1,TNFM_x)
        do j = 1,TNFM_y
          write(100,2000) j,(phi_ij(i,j,1), i=1,TNFM_x)  
        enddo
        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '             OUTPUT  PARAMETER - SOLUTION  TO  TRANSPORT  EQUATION' 
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100,3090)    'K-EFF                    =',k_eff
        write (100,3020)    'N. OUTER ITERATIONS      =',it1
        write (100,3020)    'TOTAL INNER ITERATIONS   =',it2
        write (100,4000)    'TOTAL EXECUTION TIME     =',tm,'(H:M:S)'
        write (100, FMT=* ) ''
        write (100, FMT=* ) '********************************************************************************'
        2000 format(1x,1p,i4,1x,200e13.5) 
        3000 format(1x,A8,2x,300(A14,i2))  
        3010 format(1x,A17,A22,A10)
        3020 format(1x,A26,4x,i10)
        3040 format(1x,A33,2x,200F10.5)
        3050 format(1x,1p,A41,4x,e8.1)
        3060 format(1x,1p,i11,5x,e16.5,e16.5,e16.5,e16.5)
        3070 format(1x,A18,i4)
        3080 format(1x,1p,i11,5x,e16.5,e16.5,e16.5,e16.5,e16.5)
        3090 format(1x,A26,6x,f8.6)
        4000 format(1x,A26,4x,A10,A10)
        close(100)
end subroutine Output      



