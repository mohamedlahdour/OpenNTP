
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
       pwi = pwi/4
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
             n = nfmesh_xy(j,1); m = nfmesh_xy(j,1)
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
    subroutine Flux_guess(phi_lmij,TNFM_x,TNFM_y,NORD,mu,w,order,ng)
        implicit none
        integer(kind=4), intent(in) :: TNFM_x,TNFM_y,NORD,order,ng
        real(kind=8), dimension(NORD*(NORD+2)/8), intent(in) :: mu,w
        real(kind=8), dimension(TNFM_x,TNFM_y,order,order,ng), intent(out) :: phi_lmij
        real(kind=8), dimension(TNFM_x,TNFM_y,NORD*(NORD+2)/2,ng) :: phi_ijn
        real(kind=8), dimension(TNFM_x,TNFM_y,order,order,ng) :: phi
        real(kind=8) :: plm
        integer(kind=4) :: i,j,k,l,m,n,numord,nn,mm
        numord = NORD*(NORD+2)/8
        phi = 0.0d0
        phi_ijn = 1.0d0

        do l = 0,order-1     ! itération sur l
        do m = 0,l           ! itération sur m  
        do k = 1,ng      
        do j = 1,TNFM_y
           do i = 1,TNFM_x
              nn=1
              do mm = 1,4
              do n = 1,numord
                 phi(i,j,l+1,m+1,k) = phi(i,j,l+1,m+1,k) + &
                 2.0d0*w(n)*phi_ijn(i,j,nn ,k)*plm(l,m,mu(n))
                 nn=nn+1
              enddo
              enddo
           enddo
        enddo
        enddo
        enddo
        enddo
        phi_lmij = phi
    end subroutine Flux_guess
!BL5
    subroutine Scattering_Source(TNFM_x,TNFM_y,Nmat,order,ng,NORD,&
                                 fmmid2D,mu,w,SigS,phi_lmij,Qs_ij)
        implicit none
        integer(kind=4), intent(in) :: TNFM_x,TNFM_y,Nmat,order,ng,NORD
        real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
        real(kind=8), dimension(NORD*(NORD+2)/8), intent(in) :: mu,w
        integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
        real(kind=8), dimension(TNFM_x,TNFM_y,order,order,ng), intent(in) :: phi_lmij
        real(kind=8), dimension(TNFM_x,TNFM_y,ng), intent(out) :: Qs_ij

        real(kind=8), dimension(TNFM_x,TNFM_y,order,order,ng) :: Qs_lmij
        real(kind=8), dimension(TNFM_x,TNFM_y,NORD*(NORD+2)/2,ng) :: Qs_nij
        real(kind=8), dimension(TNFM_x,TNFM_y,ng) :: Qs
        real(kind=8), parameter :: PI = 3.141592653589793
        real(kind=8) :: plm,som
        integer(kind=4) :: i,j,k,l,n,m,ii,nn,numord
        numord = NORD*(NORD+2)/8
        Qs_nij = 0.0d0; Qs = 0.0d0
        do k=1,ng
        do l = 0,order-1
           do m = 0,l
              do j = 1,TNFM_y
                 do i = 1,TNFM_x
                    Qs_lmij(i,j,l+1,m+1,k) = sum(SigS(fmmid2D(j,i),l+1,:,k)&
                                                 *phi_lmij(i,j,l+1,m+1,:))
                 enddo
              enddo
           enddo
        enddo
        enddo

        do k = 1,ng
        do j = 1,TNFM_y
           do i = 1,TNFM_x
              nn = 1
              do ii = 1,4
              do n  = 1,numord
              do l = 0,order-1
                    som = 0.0D0
                    do m = 0,l
                       som = som + Qs_lmij(i,j,l+1,m+1,k)*plm(l,m,mu(n))
                    enddo
                    Qs_nij(i,j,nn,k) = Qs_nij(i,j,nn,k)+(0.5D0)*DBLE(2*l+1)*som       
              enddo
              nn=nn+1
              enddo
              enddo
           enddo
        enddo
        enddo
        !>  update the scattering source
        do k = 1,ng
        do j = 1,TNFM_y
           do i = 1,TNFM_x
             nn=1
             do ii = 1,4
               do n  = 1,numord
               Qs(i,j,k) = Qs(i,j,k) + 2.0D0*w(n)*Qs_nij(i,j,nn,k) !> end
               nn=nn+1
               enddo
             enddo
           enddo
        enddo
        enddo
        Qs_ij = Qs
        !WRITE(*,'(11F7.3)') Qs_ij
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
        real(kind=8), dimension(TNFM_x,TNFM_y,order,order,ng), intent(in) :: phi_lmij
        real(kind=8), dimension(TNFM_x), intent(in) :: xfm
        real(kind=8), dimension(TNFM_y), intent(in) :: yfm
        real(kind=8), dimension(TNFM_x,TNFM_y,ng), intent(out) :: Qf_ij
        real(kind=8), intent(out) :: dsnew

        real(kind=8), dimension(TNFM_x,TNFM_y,NORD*(NORD+2)/2,ng) :: Qf_nij
        real(kind=8), dimension(TNFM_x,TNFM_y,order,order,ng) :: phi_11ij, Qf_lmij
        real(kind=8), dimension(TNFM_x,TNFM_y,ng) :: Qf
        real(kind=8), dimension(TNFM_x*TNFM_y*ng*order*order) :: V
        real(kind=8), parameter :: PI = 3.141592653589793                                   
        real(kind=8) :: som, plm
        integer(kind=4) :: i,j,k,l,m,n,ii,nn,n1,numord
        numord = NORD*(NORD+2)/8
        Qf_nij = 0.0d0; Qf = 0.0d0
        phi_11ij = 0.0d0
        phi_11ij(:,:,1,1,:) = phi_lmij(:,:,1,1,:)
        n1=1
        do k=1,ng
        do l = 0,order-1
           do m = 0,l
              do j = 1,TNFM_y
                 do i = 1,TNFM_x
                    Qf_lmij(i,j,l+1,m+1,k) = Chi(fmmid2D(i,j),k)*sum(NusigF(fmmid2D(i,j),:)&
                                                 *phi_11ij(i,j,1,1,:))/k_eff 
                    V(n1) = Qf_lmij(i,j,l+1,m+1,k)
                    n1 = n1+1
                 enddo
              enddo
           enddo
        enddo
        enddo
        dsnew = sum(V)

        do k = 1,ng
        do j = 1,TNFM_y
           do i = 1,TNFM_x
              nn = 1
              do ii =1,4
                 do n = 1,numord
                 do l = 0,order-1
                    som = 0.0D0
                    do m = 0,l
                     som = som + Qf_lmij(i,j,l+1,m+1,k)*plm(l,m,mu(n))
                    enddo
                    Qf_nij(i,j,nn,k) = Qf_nij(i,j,nn,k) + (0.5D0)*DBLE(2*l+1)*som
                    nn=nn+1
                 enddo
                 enddo
              enddo
           enddo
        enddo
        enddo
        !>  update the fission source

        do k = 1,ng
        do j = 1,TNFM_y
           do i = 1,TNFM_x
             nn=1
             do ii = 1,4
             do n  = 1,numord
             Qf(i,j,k) = Qf(i,j,k) + 2.0D0*w(n)*Qf_nij(i,j,nn,k) !> end
             nn=nn+1
             enddo
             enddo
           enddo
        enddo
        enddo
        Qf_ij = Qf
        print*,dsnew
        !WRITE(*,'(11F7.3)') Qf_ij
    end subroutine Fission_Source
!BL7
    subroutine Eigenvalues(TNFM_x,TNFM_y,Nmat,order,ng,N,ny,nx,&
                           fmmid2D,mu,eta,w,xfm,yfm,SigT,NusigF,&
                           BC1,BC2,BC3,BC4,Chi,SigS,phi_ij,nfmesh_xy)
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
    ! Variables Locales
    real(kind=8), dimension(TNFM_x,TNFM_y,order,order,ng) :: phi
    real(kind=8), dimension(TNFM_x,TNFM_y,N*(N+2)/2,ng) :: phi_ijn
    real(kind=8), dimension(TNFM_x,TNFM_y,order,order,ng) :: phi_lmij
    real(kind=8), dimension(TNFM_x,TNFM_y,ng) :: Qf_ij,Qs_ij,Q_ij
    real(kind=8), dimension(TNFM_x,TNFM_y,ng) :: phi_ij0,phi_ij1
    real(kind=8), dimension(TNFM_x*TNFM_y*ng) :: V
    real(kind=8) :: err_k_eff, err_phi,eps,k_eff,k_eff0,err_flux
    real(kind=8) :: dsnew, dsold,plm
    integer(kind=4) :: inter,exter,i,j,k,l,m,m2,n2,n3,numord
    call Flux_guess(phi_lmij,TNFM_x,TNFM_y,N,mu,w,order,ng)
    !call Flux_guess(phi_ij,TNFM_x,TNFM_y,NORD,mu,w,order,ng)
    phi_ij0   = phi_lmij(:,:,1,1,:)
    phi_ij1   = phi_ij
    inter     = 0
    err_k_eff = 1.0D0
    err_phi   = 1.0D0
    exter     = 0 
    k_eff     = 1.0D0 ! Given K_eff 
    k_eff0    = k_eff
    eps       = 1.0E-6
    numord    = N*(N+2)/8
    call Fission_Source(TNFM_x,TNFM_y,Nmat,order,ng,N,fmmid2D,&
                              xfm,yfm,dsnew,mu,w,NusigF,Chi,phi_lmij,k_eff,Qf_ij)

    !call Fission_Source(TNFM_x,TNFM_y,Nmat,order,ng,N,fmmid2D,&
    !                    xfm,yfm,dsnew,mu,NusigF,Chi,phi_ij,k_eff,Qf_ij)
    dsold = dsnew
    do while (err_k_eff >= eps .and. err_phi >= eps)  ! Staring extern Iteration
       call Fission_Source(TNFM_x,TNFM_y,Nmat,order,ng,N,fmmid2D,&
                              xfm,yfm,dsnew,mu,w,NusigF,Chi,phi_lmij,k_eff,Qf_ij)
       !call Fission_Source(TNFM_x,TNFM_y,Nmat,order,ng,N,fmmid2D,&
       !                       xfm,yfm,dsnew,mu,NusigF,Chi,phi_ij,k_eff,Qf_ij)
       err_flux = 1.0D0
       do while ( err_flux >= eps)  ! Starting intern Iteration
          call Scattering_Source(TNFM_x,TNFM_y,Nmat,order,ng,N,&
                                 fmmid2D,mu,w,SigS,phi_lmij,Qs_ij)
          !call Scattering_Source(TNFM_x,TNFM_y,Nmat,order,ng,N,fmmid2D,&
          !                       mu,SigS,phi_ij,Qs_ij)
          Q_ij = Qf_ij + Qs_ij
          call Sweep2D(BC1,BC2,BC3,BC4,TNFM_x,TNFM_y,N,ng,&
                       Nmat,SigT,Q_ij,mu,eta,w,xfm,yfm,fmmid2D,phi_ij,phi_ijn)
          !call Sweep2D(BC1,BC2,BC3,BC4,TNFM_x,TNFM_y,N,ng,&
          !                  Nmat,SigT,Q_ij,mu,eta,w,xfm,yfm,fmmid2D,phi_ij)    
          err_flux = maxval(abs(phi_ij-phi_ij0)/phi_ij0) ! Condition sur le flux scalaire
          phi_ij0 = phi_ij
          phi = 0.0d0
          do l = 0,order-1     ! itération sur l
          do m = 0,l           ! itération sur m  
          do k = 1,ng      
          do j = 1,TNFM_y
             do i = 1,TNFM_x
                n2=1
                do m2 = 1,4
                do n3 = 1,numord
                   phi(i,j,l+1,m+1,k) = phi(i,j,l+1,m+1,k) + &
                   2.0d0*w(n3)*phi_ijn(i,j,n2,k)*plm(l,m,mu(n3))
                   n2 = n2+1
                enddo
                enddo
             enddo
          enddo
          enddo
          enddo
          enddo
          phi_lmij = phi
          inter = inter + 1
       enddo ! ending intern Iteration    
      
       if (exter==0) go to 100    
       !print*,dsnew,dsold
       k_eff =  k_eff*dsnew/dsold
       err_k_eff =  abs(k_eff-k_eff0)/k_eff
       err_phi =  maxval(abs(phi_ij-phi_ij1)/phi_ij1)  ! Condition sur le flux scalaire
       write(*,2000)exter,k_eff,err_k_eff
       100 continue
       print*,dsnew,dsold
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
!BL9
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
    real(kind=8), dimension(TNFM_x,TNFM_y,ng), intent(out) :: phi_ij
    real(kind=8), dimension(TNFM_x,TNFM_y,NORD*(NORD+2)/2,ng), intent(out) :: phi_ijn
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
    !        mu > 0, eta > 0,   LEFT-TO-RIGHT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = 1,numord/4 
       do j = 1,TNFM_y
          phiax(1,j,n)   = 0.0D0 
          phibx(1,j,n)   = 1.0D0
          ay = 2.0D0*abs(eta(n))/dy(j)
          do i = 1,TNFM_x         
             phiay(i,1,n)   = 0.0D0
             phiby(i,1,n)   = 1.0D0
             ax = 2.0D0*abs(mu(n))/dx(i)
             deno = SigT(fmmid2D(i,j),k) + ax + ay
             nume = Q_ij(i,j,k) + ax*phiax(i,j,n) + ay*phiay(i,j,n)
             phi_aijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_aijn(i,j,n,k)<0) phi_aijn(i,j,n,k) = 0.0D0 !> end
             phiax(i+1,j,n) = 2.0D0*phi_aijn(i,j,n,k)-phiax(i,j,n)
             !> simple negative flux fixup         
             if (phiax(i+1,j,n)<0) phiax(i+1,j,n) = 0.0D0 !> end
             nume = Q_ij(i,j,k) + ax*phibx(i,j,n) + ay*phiby(i,j,n)
             phi_bijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_bijn(i,j,n,k)<0) phi_bijn(i,j,n,k) = 0.0D0 !> end
             phibx(i+1,j,n) = 2.0D0*phi_bijn(i,j,n,k)-phibx(i,j,n)
             !> simple negative flux fixup         
             if (phibx(i+1,j,n)<0) phibx(i+1,j,n) = 0.0d0 !> end
          enddo
          do i = 1,TNFM_x
             phiay(i,j+1,n) = 2.0D0*phi_aijn(i,j,n,k)-phiay(i,j,n)
             phiby(i,j+1,n) = 2.0D0*phi_bijn(i,j,n,k)-phiby(i,j,n)
             !> simple negative flux fixup         
             if (phiay(i,j+1,n)<0) phiay(i,j+1,n) = 0.0D0
             if (phiby(i,j+1,n)<0) phiby(i,j+1,n) = 0.0D0 !> end
             !> Boundary conditions in the bottom side
             if (BC3 == 'Vacuum Bottom') beta = 0.0D0
             if (BC3 == 'Reflective Bottom') beta = 1.0D0
             phiy(i,1,n)= beta*phiay(i,2,n)/(1.0D0+beta*(phiay(i,2,n)-phiby(i,2,n))) ! end
             !> simple negative flux fixup         
             if (phiy(i,1,n)<0) phiy(i,1,n) = 0.0D0 !> end
          enddo
             !> Boundary conditions in the left side
             if (BC2=='Vacuum Left') beta = 0.0D0
             if (BC2=='Reflective Left') beta = 1.0D0 
             phix(1,j,n) = beta*phiax(2,j,n)/(1.0D0+beta*(phiax(2,j,n)-phibx(2,j,n))) !> end
             !> simple negative flux fixup         
             if (phix(1,j,n)<0) phix(1,j,n) = 0.0D0 !> end
       enddo
    enddo

    do n = 1,numord/4
       do j = 1,TNFM_y
          ay = 2.0D0*abs(eta(n))/dy(j)
          do i = 1,TNFM_x
             ax = 2.0D0*abs(mu(n))/dx(i)
             deno = SigT(fmmid2D(i,j),k) + ax + ay 
             nume = Q_ij(i,j,k) + ax*phix(i,j,n) + ay*phiy(i,j,n)
             phi_ijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0D0 !> end
             phix(i+1,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i,j,n)
             !> simple negative flux fixup         
             if (phix(i+1,j,n)<0) phix(i+1,j,n) = 0.0D0 !> end
             !>  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 2.0D0*w(n)*phi_ijn(i,j,n,k) !> end
          enddo
          do i = 1,TNFM_x
             phiy(i,j+1,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j,n)
             !> simple negative flux fixup         
             if (phiy(i,j+1,n)<0) phiy(i,j+1,n) = 0.0D0 !> end
          enddo
       enddo
    enddo
    enddo

    !       -----------------------------------------------------------------        
    !        mu < 0, eta < 0,   RIGHT-TO-LEFT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    do k = 1,ng
    nn=1
    do n = numord/4+1,numord/2
       do j = TNFM_y,1,-1
          phiax(TNFM_x+1,j,n) = 0.0D0 
          phibx(TNFM_x+1,j,n) = 1.0D0
          ay = 2.0d0*abs(eta(nn))/dy(j)
          do i = TNFM_x,1,-1         
             phiay(i,TNFM_y+1,n) = 0.0D0
             phiby(i,TNFM_y+1,n) = 1.0D0
             ax = 2.0d0*abs(mu(nn))/dx(i)
             deno = SigT(fmmid2D(i,j),k) + ax + ay
             nume = Q_ij(i,j,k) + ax*phiax(i+1,j,n) + ay*phiay(i,j+1,n)
             phi_aijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_aijn(i,j,n,k)<0) phi_aijn(i,j,n,k) = 0.0d0 !> end
             phiax(i,j,n) = 2.0D0*phi_aijn(i,j,n,k)-phiax(i+1,j,n)
             !> simple negative flux fixup         
             if (phiax(i,j,n)<0) phiax(i,j,n) = 0.0D0 !> end
             nume = Q_ij(i,j,k) + ax*phibx(i+1,j,n) + ay*phiby(i,j+1,n)
             phi_bijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_bijn(i,j,n,k)<0) phi_bijn(i,j,n,k) = 0.0D0 !> end 
             phibx(i,j,n) = 2.0D0*phi_bijn(i,j,n,k)-phibx(i+1,j,n)
             !> simple negative flux fixup         
             if (phibx(i,j,n)<0) phibx(i,j,n) = 0.0D0 !> end
          enddo
          do i = TNFM_x,1,-1         
             phiay(i,j,n) = 2.0D0*phi_aijn(i,j,n,k)-phiay(i,j+1,n)
             phiby(i,j,n) = 2.0D0*phi_bijn(i,j,n,k)-phiby(i,j+1,n)
             !> simple negative flux fixup         
             if (phiay(i,j,n)<0) phiay(i,j,n) = 0.0d0
             if (phiby(i,j,n)<0) phiby(i,j,n) = 0.0d0 !> end
             if (BC4 == 'Vacuum Top') beta = 0.0D0
             if (BC4 == 'Reflective Top') beta = 1.0D0
             phiy(i,TNFM_y+1,n) = beta*phiay(i,TNFM_y,n)/&
                 (1.0D0+beta*(phiay(i,TNFM_y,n)-phiby(i,TNFM_y,n))) !> end
             !> simple negative flux fixup         
             if (phiy(i,TNFM_y+1,n)<0) phiy(i,TNFM_y+1,n) = 0.0D0 !> end
          enddo
             !> Boundary conditions in the right side
             if (BC1=='Vacuum Right') beta = 0.0D0
             if (BC1=='Reflective Right') beta = 1.0D0
             phix(TNFM_x+1,j,n) = beta*phiax(TNFM_x,j,n)/&
                 (1.0D0+beta*(phiax(TNFM_x,j,n)-phibx(TNFM_x,j,n))) ! end
             !> simple negative flux fixup         
             if (phix(TNFM_x+1,j,n)<0) phix(TNFM_x+1,j,n) = 0.0D0 !> end
       enddo
       nn=nn+1
    enddo   
    nn=1
    do n = numord/4+1,numord/2 
       do j = TNFM_y,1,-1  
          ay = 2.0d0*abs(eta(nn))/dy(j)       
          do i = TNFM_x,1,-1         
             ax = 2.0d0*abs(mu(nn))/dx(i)
             deno = SigT(fmmid2D(i,j),k) + ax + ay 
             nume = Q_ij(i,j,k) + ax*phix(i+1,j,n) + ay*phiy(i,j+1,n)
             phi_ijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0 !> end
             phix(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i+1,j,n)
             !> simple negative flux fixup         
             if (phix(i,j,n)<0) phix(i,j,n) = 0.0d0 !> end
             !>  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 2.0D0*w(nn)*phi_ijn(i,j,n,k) !> end
          enddo
          do i = TNFM_x,1,-1         
             phiy(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j+1,n)
             !> simple negative flux fixup         
             if (phiy(i,j,n)<0) phiy(i,j,n) = 0.0d0 !> end
          enddo
       enddo
       nn=nn+1
    enddo
    enddo
    !       -----------------------------------------------------------------        
    !        mu > 0, eta < 0,   LEFT-TO-RIGHT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    do k = 1,ng
    nn=1
    do n = numord/2+1,3*numord/4
       do j = TNFM_y,1,-1
          phiax(1,j,n) = 0.0D0 
          phibx(1,j,n) = 1.0D0
          ay = 2.0D0*abs(eta(nn))/dy(j)
          do i = 1,TNFM_x         
             phiay(i,TNFM_y+1,n) = 0.0D0
             phiby(i,TNFM_y+1,n) = 1.0D0
             ax = 2.0D0*abs(mu(nn))/dx(i)
             deno = SigT(fmmid2D(i,j),k) + ax + ay
             nume = Q_ij(i,j,k) + ax*phiax(i,j,n) + ay*phiay(i,j+1,n)
             phi_aijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_aijn(i,j,n,k)<0) phi_aijn(i,j,n,k) = 0.0D0 !> end
             phiax(i+1,j,n) = 2.0D0*phi_aijn(i,j,n,k)-phiax(i,j,n)
             !> simple negative flux fixup         
             if (phiax(i+1,j,n)<0) phiax(i+1,j,n) = 0.0D0 !> end
             nume = Q_ij(i,j,k) + ax*phibx(i,j,n) + ay*phiby(i,j+1,n)
             phi_bijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_bijn(i,j,n,k)<0) phi_bijn(i,j,n,k) = 0.0D0 !> end
             phibx(i+1,j,n) = 2.0D0*phi_bijn(i,j,n,k)-phibx(i,j,n)
             !> simple negative flux fixup         
             if (phibx(i+1,j,n)<0) phibx(i+1,j,n) = 0.0D0 !> end
          enddo
          do i = 1,TNFM_x
             phiay(i,j,n) = 2.0D0*phi_aijn(i,j,n,k)-phiay(i,j+1,n)
             phiby(i,j,n) = 2.0D0*phi_bijn(i,j,n,k)-phiby(i,j+1,n)
             !> simple negative flux fixup         
             if (phiay(i,j,n)<0) phiay(i,j,n) = 0.0D0
             if (phiby(i,j,n)<0) phiby(i,j,n) = 0.0D0
             !> Boundary condition in the top side
             if (BC4=='Vacuum Top') beta = 0.0D0
             if (BC4=='Reflective Top') beta = 1.0D0
             phiy(i,TNFM_y+1,n) = beta*phiay(i,TNFM_y,n)/&
                 (1.0D0+beta*(phiay(i,TNFM_y,n)-phiby(i,TNFM_y,n))) !> end
             !> simple negative flux fixup         
             if (phiy(i,TNFM_y+1,n)<0) phiy(i,TNFM_y+1,n) = 0.0D0 !> end
          enddo
             !> Boundary condition in the left side
             if (BC2=='Vacuum Left') beta = 0.0D0
             if (BC2=='Reflective Left') beta = 1.0D0
             phix(1,j,n) = beta*phiax(2,j,n)/(1.0D0+beta*(phiax(2,j,n)-phibx(2,j,n))) !> end
             !> simple negative flux fixup         
             if (phix(1,j,n)<0) phix(1,j,n) = 0.0D0 !> end
       enddo
       nn=nn+1
    enddo

    nn=1
    do n = numord/2+1,3*numord/4
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
             phi(i,j,k) = phi(i,j,k) + 2.0D0*w(nn)*phi_ijn(i,j,n,k) !> end
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
    do k = 1,ng
    nn=1
    do n = 3*numord/4+1,numord
       do j = 1,TNFM_y
          phiax(TNFM_x+1,j,n) = 0.0D0 
          phibx(TNFM_x+1,j,n) = 1.0D0
          ay = 2.0d0*abs(eta(nn))/dy(j)
          do i = TNFM_x,1,-1  
             phiay(i,1,n) = 0.0D0
             phiby(i,1,n) = 1.0D0     
             ax = 2.0d0*abs(mu(nn))/dx(i)
             deno = SigT(fmmid2D(i,j),k) + ax + ay
             nume = Q_ij(i,j,k) + ax*phiax(i+1,j,n) + ay*phiay(i,j,n)
             phi_aijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_aijn(i,j,n,k)<0) phi_aijn(i,j,n,k) = 0.0D0 !> end
             phiax(i,j,n) = 2.0D0*phi_aijn(i,j,n,k)-phiax(i+1,j,n)
             !> simple negative flux fixup         
             if (phiax(i,j,n)<0) phiax(i,j,n) = 0.0D0 !> end
             nume = Q_ij(i,j,k) + ax*phibx(i+1,j,n) + ay*phiby(i,j,n)
             phi_bijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_bijn(i,j,n,k)<0) phi_bijn(i,j,n,k) = 0.0D0 !> end
             phibx(i,j,n) = 2.0D0*phi_bijn(i,j,n,k)-phibx(i+1,j,n)
             !> simple negative flux fixup         
             if (phibx(i,j,n)<0) phibx(i,j,n) = 0.0D0 !> end
          enddo
          do i = TNFM_x,1,-1         
             phiay(i,j+1,n) = 2.0D0*phi_aijn(i,j,n,k)-phiay(i,j,n)
             phiby(i,j+1,n) = 2.0D0*phi_bijn(i,j,n,k)-phiby(i,j,n)
             !> simple negative flux fixup         
             if (phiay(i,j+1,n)<0) phiay(i,j+1,n) = 0.0D0
             if (phiby(i,j+1,n)<0) phiby(i,j+1,n) = 0.0D0 !> end
             !> Boundary condition in the bottom side
             if (BC3 == 'Vacuum Bottom') beta = 0.0D0
             if (BC3 == 'Reflective Bottom') beta = 1.0D0
             phiy(i,1,n) = beta*phiay(i,2,n)/(1.0D0+beta*(phiay(i,2,n)-phiby(i,2,n)))
             !> simple negative flux fixup         
             if (phiy(i,1,n)<0) phiy(i,1,n) = 0.0D0 !> end
          enddo
             !> Boundary condition in the right side
             if (BC1 == 'Vacuum Right') beta = 0.0D0
             if (BC1 == 'Reflective Right') beta = 1.0D0
             phix(TNFM_x+1,j,n) = beta*phiax(TNFM_x,j,n)/&
                 (1.0D0+beta*(phiax(TNFM_x,j,n)-phibx(TNFM_x,j,n))) !> end
             !> simple negative flux fixup         
             if (phix(TNFM_x+1,j,n)<0) phix(TNFM_x+1,j,n) = 0.0D0 !> end
       enddo
       nn=nn+1
    enddo
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
             phi(i,j,k) = phi(i,j,k) + 2.0D0*w(nn)*phi_ijn(i,j,n,k)
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
    phi_ij = phi
    end subroutine Sweep2D
!BL10
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
       print*, '_____________________________________________________________'
       print*, '          | The OpenRSN Neutron Transport Package            '       
       print*, 'Version   | Version Number: 1.2                              '     
       print*, 'Copyright | 2015-2019 Radiation and Nuclear system Laboratory'
       print*, '          | University Abdelmalk Essaadi, FS Tetouan, Morocco'
       print*, 'Source    | FORTRAN90 version                                ' 
       print*, 'GUI       | PyQt5                                            ' 
       print*, 'Method    | The Discrete Ordinates Method Sn                 '  
       print*, 'Dimension | Two dimensions (2D)                              '
       print*, 'Geometry  | Cartesian                                        ' 
       print*, '_____________________________________________________________'
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


