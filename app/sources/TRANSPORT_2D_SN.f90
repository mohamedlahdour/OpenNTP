!>BLOC1
    subroutine quadrature_set(N,mup, etap, psip,pwi)
       !> Level Symmetric Gaussian Quadrature Sets
       !> The level-symmetric quadrature set is used in the Discrete Ordinates (SN) 
       !> method (Lewis & Miller).The subscript refers to the number of directions 
       !> along each !axis with half being positive and half negative. 
       !> This subroutine below give the weights and angles used for each set in the 1st octant.
       implicit none
       integer(kind=4), intent(in) :: N
       real(kind=4), dimension(N/2) :: mmu
       real(kind=4), dimension(N*(N+2)/8), intent(out) :: mup, etap, psip,pwi
       real(kind=4), parameter :: pi = 3.141592653589793
       real(kind=4) :: mu1,C
       integer(kind=4) :: i,j,k,ip,M
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
           C = 2.0*(1.0-3.0*mu1)/DBLE(N-2)
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
          ! if k is not 0 it means that line contains "CR" 
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
!>BLOC2
    subroutine fmm_id2D(core,assembly,fmmid2D,nfmesh_xy,RegMat,npx,npy,npc,tnfm_x,tnfm_y,nx,ny,nyy,nxx,na)
       !> This subroutine fmm_id2D(i,j) allowing to identify which type of material in the region (i, j)
       ! -- fine mesh material id x and y directions 
       !>  npc        Total number of PinCell
       !>  nxx        Total number of Coarse mesh along x axis for each assembly.
       !>  nyy        Total number of Coarse mesh along y axis for each assembly.
       !>  na         Total number of assembly
       !>  nx         Total number of Coarse mesh along x axis.
       !>  ny         Total number of Coarse mesh along y axis.
       !>  npx        Total number of Coarse mesh along x axis for each pincell.
       !>  npy        Total number of Coarse mesh along y axis for each pincell.
       !>  TNFM_x     Total number of fine mesh along x axis.
       !>  TNFM_y     Total number of fine mesh along y axis.
       !>  RegMat     Which material in pincell with index i (i = 1, ..., npc) along npx axis and npy axis
       !>  nfmesh_xy  Numbre of finite mesh in each pin cell along npx axis and npy axis
       implicit none
       integer(kind=4), intent(in) :: npx,npy,npc,tnfm_x,tnfm_y,nx,ny,nyy,nxx,na
       integer(kind=4), dimension(nx,ny), intent(in) :: core
       integer(kind=4), dimension(na,nxx,nyy), intent(in) :: assembly
       integer(kind=4), dimension(npc,npx,npy), intent(in) :: RegMat, nfmesh_xy
       integer(kind=4), dimension(tnfm_x,tnfm_y), intent(out) :: fmmid2D
       ! Variables Locales
       integer(kind=4) :: n1,n2,n3,n4,n5,n6,m1,m2,somx,somy
       m1=1
       do n1=1,ny
       do n3=1,nyy
       do n6 = 1,npy
       m2=1
       do n2 = 1,nx
          do n4 = 1,nxx
              do n5 = 1,npx      
                somx = nfmesh_xy(assembly(core(n2,n1),n4,n3),n5,n6)
                somy = nfmesh_xy(assembly(core(n2,n1),n4,n3),n5,n6)
                fmmid2D(m2:m2+somx-1,m1:m1+somy-1) = RegMat(assembly(core(n2,n1),n4,n3),n5,n6)
              m2=m2+somx
              enddo
          enddo
       enddo
       m1=m1+somy
       enddo
       enddo
       enddo
    end subroutine fmm_id2D
!>BLOC3
    subroutine Mesh_2D(nx,ny,npc,npx,npy,nxx,nyy,na,TNFM_x,TNFM_y,core,assembly,nfmesh_xy,xcm,ycm,xfm,yfm)
       !>  npc       Total number of PinCell
       !>  nx        Total number of Coarse mesh along x axis.
       !>  ny        Total number of Coarse mesh along y axis.
       !>  npx       Total number of Coarse mesh along x axis for each pincell.
       !>  npy       Total number of Coarse mesh along y axis for each pincell.
       !>  TNFM_x    Total number of fine mesh along x axis.
       !>  TNFM_y    Total number of fine mesh along y axis.
       !>  xfm       Number of fine mesh per x coarse mesh for each pincell.
       !>  yfm       Number of fine mesh per y coarse mesh for each pincell.
       !>  xcm       Coarse mesh edges along x axis.
       !>  ycm       Coarse mesh edges along y axis.
       !>  Size each mesh in the x direction
       implicit none
       integer(kind=4), intent(in) :: npx,npy,TNFM_x,TNFM_y,nx,ny,npc,na,nxx,nyy
       integer(kind=4), dimension(npc,npx,npy), intent(in) ::  nfmesh_xy
       real(kind=4), dimension(npc,npx), intent(in)  :: xcm
       real(kind=4), dimension(npc,npy), intent(in)  :: ycm
       integer(kind=4), dimension(nx,ny), intent(in) :: core
       integer(kind=4), dimension(na,nxx,nyy), intent(in) :: assembly
       real(kind=4), dimension(TNFM_x), intent(out)  :: xfm
       real(kind=4), dimension(TNFM_y), intent(out)  :: yfm
       integer(kind=4) :: n1,n2,n3,n,som
       !>  Size each mesh in the x direction
       n=1
       do n1 = 1,nx
          do n2 = 1,nxx
              do n3 = 1,npx 
                som = nfmesh_xy(assembly(core(n1,1),n2,1),n3,1)
                xfm(n:n+som-1) = xcm(assembly(core(n1,1),n2,1),n3)/nfmesh_xy(assembly(core(n1,1),n2,1),n3,1)
                n=n+som
              enddo
          enddo
       enddo
      !>  Size each mesh in the y direction
       n=1
       do n1 = 1,ny
          do n2 = 1,nyy
              do n3 = 1,npy 
                som = nfmesh_xy(assembly(core(1,n1),1,n2),1,n3)
                yfm(n:n+som-1) = ycm(assembly(core(1,n1),1,n2),n3)/nfmesh_xy(assembly(core(1,n1),1,n2),1,n3)
                n=n+som
              enddo
          enddo
       enddo
    end subroutine Mesh_2D
!>BLOC4
    subroutine GuessFlux(phi_ij,TNFM_x,TNFM_y,xfm,yfm,NusigF,Chi,fmmid2D,&
                         order,Nmat,ng,N)
       !> This function allows to calculate the guess flux
       !> phi_ij    : The guess SCALAR flux
       !> Delta     : The volume of each cell
       !> NusigF    : NuFission cross section
       !> xcm       : Coarse mesh edges along x axis.
       !> ycm       : Coarse mesh edges along y axis.
        implicit none
        integer(kind=4), intent(in) :: TNFM_x,TNFM_y,order,Nmat,ng,N
        integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
        real(kind=4), dimension(Nmat,ng), intent(in) :: NusigF,Chi
        real(kind=4), dimension(TNFM_x), intent(in) :: xfm
        real(kind=4), dimension(TNFM_y), intent(in) :: yfm
        real(kind=4), dimension(TNFM_x,TNFM_y,ng), intent(out) :: phi_ij
        real(kind=4) :: val
        integer(kind=4) :: i,j,k,k1     
           do j = 1,TNFM_y
              do i = 1,TNFM_x
                 val = 0.
                 do k1 = 1,Nmat
                    val = val + xfm(i)*yfm(j)*sum(NusigF(k1,:))
                 enddo
                 phi_ij(i,j,:) = 1.0/val
              enddo
           enddo  
    end subroutine GuessFlux
!>BLOC5
    subroutine Scattering_Source(TNFM_x,TNFM_y,Nmat,order,ng,N,fmmid2D,SigS,phi_ij,Qs_ijn)
       !> Calculate the scatternig source components Qs_ijn.
       !> TNFM_x  : Total number of fine mesh along x axis.
       !> TNFM_y  : Total number of fine mesh along y axis.
       !> ng      : The total number of energy group
       !> SigS    : Scatternig cross section
       !> phi_ij  : The scalar flux
       !> N       : Number of Discrete Ordinates
        implicit none
        integer(kind=4), intent(in) :: TNFM_x,TNFM_y,Nmat,order,ng,N
        real(kind=4), dimension(Nmat,order,ng,ng), intent(in) :: SigS
        integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
        real(kind=4), dimension(TNFM_x,TNFM_y,ng), intent(in) :: phi_ij
        real(kind=4), dimension(TNFM_x,TNFM_y,N*(N+2)/2,ng), intent(out) :: Qs_ijn
        ! Local Variables 
        integer(kind=4) :: i,j,k
        do k = 1,ng
        do j = 1,TNFM_y
        do i = 1,TNFM_x
           Qs_ijn(i,j,:,k) = 0.25*sum(SigS(fmmid2D(i,j),1,:,k)*phi_ij(i,j,:))
        enddo
        enddo
        enddo
    end subroutine Scattering_Source
!>BLOC6
    subroutine Fission_Source(TNFM_x,TNFM_y,Nmat,order,ng,N,fmmid2D,dsnew,NusigF,Chi,phi_ij,k_eff,Qf_ijn)
        !> Calculate the Fission source components Qf_ijn.
        !> TNFM_x  : Total number of fine mesh along x axis.
        !> TNFM_y  : Total number of fine mesh along y axis.
        !> ng      : The total number of energy group
        !> NuSigF  : NuFission cross section
        !> CHI     : Density function for neutrons
        !> phi_ij  : The scalar flux
        !> N       : Number of Discrete Ordinates
        !> k_eff   : Multiplication factor
        implicit none
        integer(kind=4), intent(in) :: TNFM_x,TNFM_y,Nmat,order,ng,N
        real(kind=4), intent(in) :: k_eff 
        real(kind=4), dimension(Nmat,ng), intent(in) :: NusigF,Chi
        integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
        real(kind=4), dimension(TNFM_x,TNFM_y,ng), intent(in) :: phi_ij
        real(kind=4), dimension(TNFM_x,TNFM_y,N*(N+2)/2,ng), intent(out) :: Qf_ijn
        real(kind=4), intent(out) :: dsnew
        ! Varaibles locales
        integer(kind=4) :: i,j,k
        !>--------------- 
        do k = 1,ng
        do j = 1,TNFM_y
        do i = 1,TNFM_x
           Qf_ijn(i,j,:,k) = 0.25*Chi(fmmid2D(i,j),k)*&
              sum(NusigF(fmmid2D(i,j),:)*phi_ij(i,j,:))/k_eff
        enddo
        enddo
        enddo   
        dsnew = sum(Qf_ijn)
        !>--------------- 
    end subroutine Fission_Source
!>BLOC7
    subroutine Eigenvalues(TNFM_x,TNFM_y,Nmat,order,ng,N,ny,nx,eps,Max_it,&
                           fmmid2D,mu,eta,psi,w,xfm,yfm,SigT,NusigF,SigF,&
                           BC1,BC2,BC3,BC4,Chi,SigS,phi_ij,nfmesh_xy,k_eff,inter,exter)
    !> Iterative method to calculate eigenvalues and vectors
    !> phi_ij     : Scalar flux
    !> Max_it     : Maximum number of iterations
    !> k_eff      : multiplication factor
    !> eps        : Criterion of convergence
    !> inter      : Number of internal iterations
    !> exter      : Number of external iterations
    !> SigF       : Fission cross section
    !> NuSigF     : NuFission cross section
    !> CHI        : Density function for neutrons
    !> SigS       : Scatternig cross section
    !> SigT       : Total cross section
    !> mu,eta,psi : the direction cosines
    implicit none
    ! Global Variables ..........................................
    CHARACTER(50), intent(in) :: BC1,BC2,BC3,BC4
    real(kind=4), intent(in) :: eps
    integer(kind=4), intent(in) :: TNFM_x,TNFM_y,Nmat,order,ng,N,ny,nx,Max_it
    integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
    real(kind=4), dimension(N*(N+2)/8), intent(in) :: mu,eta,psi,w
    real(kind=4), dimension(TNFM_x), intent(in) :: xfm
    real(kind=4), dimension(TNFM_y), intent(in) :: yfm
    real(kind=4), dimension(Nmat,ng), intent(in) :: SigT,NusigF,Chi,SigF
    real(kind=4), dimension(Nmat,order,ng,ng), intent(in) :: SigS
    integer(kind=4), dimension(ny,nx), intent(in) :: nfmesh_xy
    real(kind=4), dimension(TNFM_x,TNFM_y,ng), intent(out) :: phi_ij
    real(kind=4), intent(out) :: k_eff
    integer(kind=4), intent(out) :: inter,exter
    ! Local Variables ..........................................
    real(kind=4), dimension(TNFM_x,TNFM_y,N*(N+2)/2,ng) :: phi_ijn
    real(kind=4), dimension(TNFM_x,TNFM_y,N*(N+2)/2,ng) :: Qf_ijn,Qs_ijn,Q_ijn
    real(kind=4), dimension(TNFM_x,TNFM_y,ng) :: phi_ij0,phi_ij1
    real(kind=4), dimension(TNFM_x*TNFM_y*ng) :: V,V1
    real(kind=4) :: err_k_eff, err_phi,k_eff0,err_flux
    real(kind=4) :: dsnew, dsold,Rlm,err1,err2
    integer(kind=4) :: i,j,k,l,m,n1,i1,n2,mm
    !> -------------------------------------
    call GuessFlux(phi_ij,TNFM_x,TNFM_y,xfm,yfm,NusigF,Chi,fmmid2D,order,Nmat,ng,N) 
    err_k_eff = 1.0
    err_phi   = 1.0
    exter     = 1
    k_eff     = 1.0 ! Given K_eff 
    k_eff0    = 1.0
    call Fission_Source(TNFM_x,TNFM_y,Nmat,order,ng,N,fmmid2D,&
                              dsnew,NusigF,Chi,phi_ij,k_eff,Qf_ijn)
    dsold = dsnew
    do while (err_k_eff >= eps .and. err_phi >= eps)  ! Staring extern Iteration
       phi_ij1 = phi_ij
       err_flux = 1.0
           if ( exter > Max_it ) then
           print*, 'error( unable to converge )'
           stop
           end if
           inter = 1
       do while (err_flux >= eps*10)  ! Starting intern Iteration
          phi_ij0 = phi_ij
          call Scattering_Source(TNFM_x,TNFM_y,Nmat,order,ng,N,fmmid2D,SigS,phi_ij,Qs_ijn)
          Q_ijn = Qf_ijn + Qs_ijn 
          call Sweep2D(Nmat,ng,TNFM_x,TNFM_y,N,Q_ijn,xfm,yfm,SigT,&
                     mu,eta,w,fmmid2D,phi_ij)
          ! Condition on scalar flux during internal iteration
          err1   = maxval(abs(phi_ij0))
          err2   = maxval(abs(phi_ij-phi_ij0))
          err_flux = err2/err1
          if (inter == 5000) exit
          inter = inter + 1     
       enddo ! ending intern Iteration
       ! Condition on scalar flux during external iteration
       err1   = maxval(abs(phi_ij1))
       err2   = maxval(abs(phi_ij-phi_ij1))
       err_phi =  err2/err1
       call Fission_Source(TNFM_x,TNFM_y,Nmat,order,ng,N,fmmid2D,&
                              dsnew,NusigF,Chi,phi_ij,k_eff,Qf_ijn)
       ! Normalized flux
       call NormalizeFlux(TNFM_x,TNFM_y,Nmat,ng,sigF,Chi,fmmid2D,xfm,yfm,phi_ij)
       k_eff =  k_eff*dsnew/dsold
       err_k_eff =  abs(k_eff-k_eff0)/k_eff0
       write(*,2000)exter,inter,k_eff,err_k_eff
       dsold = dsnew
       k_eff0   = k_eff
       exter    = exter + 1
    enddo ! ending extern Iteration  
    2000 format("ExtIter",i4,":",1x,"IntIter",i5,":"1x,"keff =",F9.6,1x,"res =",e10.3)
    end subroutine Eigenvalues
!>BLOC8
    subroutine Plot(TNFM_x,TNFM_y,ny,nx,npx,npy,npc,ng,Nmat,napc,na,nxx,nyy,xfm,yfm,&
                    core,assembly,NusigF,SigF,Chi,fmmid2D,nfmesh_xy,phi_ij)
       !> Plot scalar flux
       implicit none
       integer(kind=4), intent(in) :: TNFM_x,TNFM_y,ny,nx,npx,npy,npc,ng,Nmat,napc,na,nxx,nyy
       real(kind=4), dimension(TNFM_x), intent(in) :: xfm
       real(kind=4), dimension(TNFM_y), intent(in) :: yfm
       real(kind=4), dimension(Nmat,ng), intent(in) :: NusigF,Chi,SigF
       integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
       real(kind=4), dimension(TNFM_x,TNFM_y,ng), intent(in) :: phi_ij
       integer(kind=4), dimension(npc,npx,npy), intent(in) ::  nfmesh_xy
       integer(kind=4), dimension(nx,ny), intent(in) :: core
       integer(kind=4), dimension(na,nxx,nyy), intent(in) :: assembly
       real(kind=4), dimension(nx*nxx,ny*nyy) :: P,PD
       real(kind=4), dimension(TNFM_y) :: dy
       real(kind=4), dimension(TNFM_x) :: dx
       real(kind=4) :: r
       integer(kind=4) :: i,j,n,k,m,k1,som1,som2,n1,n2,n3,n4,som
       ! The part common to all files
       character(14), dimension(11:10+ng) :: Filename ! The name of the complete file
       open (1000,file='app/Output/PinPower.h')
       open (1001,file='app/Output/PowerDistribution.h')

       do i=11,10+ng ! loop through all files
       write(Filename(i),'(a13,i1)') 'app/Output/FG',i-10
       enddo

       m=1;r=0.0
       do n1 = 1,ny
          do n2 = 1,nyy
              do n3 = 1,npy 
                do n4 = 1,nfmesh_xy(assembly(core(1,n1),1,n2),1,n3)
                    dy(m) =sum(yfm(1:m))
                    m=m+1
                enddo
              enddo
          enddo
       enddo

       m=1
       do n1 = 1,nx
          do n2 = 1,nxx
              do n3 = 1,npx 
                do n4 = 1,nfmesh_xy(assembly(core(n1,1),n2,1),n3,1)
                   dx(m) = sum(xfm(1:m))
                   m=m+1
                enddo
              enddo
          enddo
       enddo

       call PinPower(TNFM_y,TNFM_x,Nmat,ng,nx,ny,npx,npy,npc,napc,na,nxx,nyy,xfm,yfm,&
                            core,assembly,nfmesh_xy,phi_ij,fmmid2D,sigF,P,PD)
       som1  = sum(nfmesh_xy(assembly(core(1,1),1,1),:,1)) 
       som2  = sum(nfmesh_xy(assembly(core(1,1),1,1),1,:)) 
       write(1000,2000) r,(sum(xfm(1:som1*i)), i=1,nx*nxx) 
       write(1001,2000) r,(sum(xfm(1:som1*i)), i=1,nx*nxx) 
       do j = 1,nyy*ny
          write(1000,2000) sum(yfm(1:som2*j)),(P(i,j), i=nx*nxx,1,-1)  
       enddo
       do j = 1,nyy*ny
          write(1001,2000) sum(yfm(1:som2*j)),(PD(i,j), i=nx*nxx,1,-1)  
       enddo
       close(1000)
       close(1001)

       do k=11,10+ng
          open (k,file=Filename(k))
          write(k,2000) r,(dx(i), i=1,TNFM_x)
          do j = 1,TNFM_y
             write(k,2000) dy(j),(phi_ij(i,j,k-10), i=TNFM_x,1,-1)
          enddo
          close(k)
       enddo
       2000 format(f15.8,10000f15.8)
    end subroutine Plot
!>BLOC9
    subroutine Sweep2D(Nmat,ng,TNFM_x,TNFM_y,NORD,Q_ijn,dx,dy,SigT,&
                     mu,eta,w,fmmid2D,phi_ij)
    !> Sweep over a 2D mesh.
    !> Q_ijn   : Total source components
    !> phi_ij  : The scalar flux
    !> dx      : Number of fine mesh per x coarse mesh for each pincell.
    !> dy      : Number of fine mesh per y coarse mesh for each pincell.
    !> NORD    : Number of Discrete Ordinates
    !> psiV    : The vertical flux components
    !> psiH    : The horizontal flux components
    !> mu, eta : The direction cosines
    !> w       : The weights
    implicit none
    ! Varaibles globales
    integer(kind=4), intent(in) :: Nmat,ng,TNFM_x,TNFM_y,NORD
    integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
    real(kind=4), dimension(TNFM_x,TNFM_y,NORD*(NORD+2)/2,ng), intent(in) :: Q_ijn 
    real(kind=4), dimension(NORD*(NORD+2)/8), intent(in) :: mu, eta, w
    real(kind=4), dimension(TNFM_x), intent(in) :: dx
    real(kind=4), dimension(TNFM_y), intent(in) :: dy
    real(kind=4), dimension(Nmat,ng), intent(in) :: SigT
    real(kind=4), dimension(TNFM_x,TNFM_y,ng), intent(out) :: phi_ij
    ! Variables locales
    real(kind=4), dimension(TNFM_x+1,TNFM_y,NORD*(NORD+2)/2,ng) :: psiH
    real(kind=4), dimension(TNFM_x,TNFM_y+1,NORD*(NORD+2)/2,ng) :: psiV
    real(kind=4), dimension(TNFM_x,TNFM_y,ng) :: phi
    real(kind=4) :: ax,ay,nume,deno,psiC
    integer(kind=4) :: i,j,n,k,numord,nn
    psiH=0.0;psiV=0.0;psiC=0.0;phi=0.0
    numord = NORD*(NORD+2)/2
    !       -----------------------------------------------------------------        
    !        mu < 0, eta < 0,   RIGHT-TO-LEFT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    psiH(TNFM_x+1,TNFM_y:1:-1,numord/4:1:-1,:) = 0.0
    psiV(TNFM_x:1:-1,TNFM_y+1,numord/4:1:-1,:) = 0.0 
    do k = 1,ng
    do n = 1,numord/4
       do j = TNFM_y,1,-1
             ay = 2.0*eta(n)/dy(j)
          do i = TNFM_x,1,-1
             ax = 2.0*mu(n)/dx(i)
             nume = Q_ijn(i,j,n,k) + ax*psiH(i+1,j,n,k) + ay*psiV(i,j+1,n,k)
             deno = ax + ay + SigT(fmmid2D(i,j),k)
             psiC = nume/deno
             psiH(i,j,n,k) = 2.0*psiC-psiH(i+1,j,n,k)
             psiV(i,j,n,k) = 2.0*psiC-psiV(i,j+1,n,k)
             !  Negative flux fix up
             if (psiH(i,j,n,k)<0.0) then
                 psiH(i,j,n,k) = 0.0
             endif
             if (psiV(i,j,n,k)<0.0) then
                 psiV(i,j,n,k) = 0.0
             endif
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + w(n)*psiC
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------        
    !        mu > 0, eta < 0,   LEFT-TO-RIGHT, TOP-TO-BOTTOM
    !       ----------------------------------------------------------------- 
    psiH(1,TNFM_y:1:-1,numord/4+1:numord/2,:) = psiH(1,TNFM_y:1:-1,1:numord/4,:) 
    psiV(1:TNFM_x,TNFM_y+1,numord/4+1:numord/2,:) = 0.0 
    do k = 1,ng
    nn=1
    do n = numord/4+1,numord/2
       do j = TNFM_y,1,-1
             ay = 2.0*eta(nn)/dy(j)
          do i = 1,TNFM_x
             ax = 2.0*mu(nn)/dx(i)
             nume = Q_ijn(i,j,n,k) + ax*psiH(i,j,n,k) + ay*psiV(i,j+1,n,k)
             deno = ax + ay + SigT(fmmid2D(i,j),k)
             psiC = nume/deno
             psiH(i+1,j,n,k) = 2.0*psiC-psiH(i,j,n,k)
             psiV(i,j,n,k) = 2.0*psiC-psiV(i,j+1,n,k)
             !  Negative flux fix up
             if (psiH(i+1,j,n,k)<0) then
                 psiH(i+1,j,n,k) = 0
             endif
             if (psiV(i,j,n,k)<0) then
                 psiV(i,j,n,k) = 0
             endif
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + w(nn)*psiC
          enddo
       enddo
       nn=nn+1
    enddo
    enddo
    !       -----------------------------------------------------------------
    !        mu < 0, eta > 0,   RIGHT-TO-LEFT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    psiH(TNFM_x+1,1:TNFM_y,numord/2+1:3*numord/4,:)  = 0.0 
    psiV(TNFM_x:1:-1,1,numord/2+1:3*numord/4,:) = psiV(TNFM_x:1:-1,1,1:numord/4,:) 
    do k = 1,ng
    nn=1
    do n = numord/2+1,3*numord/4 
       do j = 1,TNFM_y
             ay = 2.0*eta(nn)/dy(j)
          do i = TNFM_x,1,-1
             ax = 2.0*mu(nn)/dx(i)
             nume = Q_ijn(i,j,n,k) + ax*psiH(i+1,j,n,k) + ay*psiV(i,j,n,k)
             deno = ax + ay + SigT(fmmid2D(i,j),k)
             psiC = nume/deno
             psiH(i,j,n,k) = 2.0*psiC-psiH(i+1,j,n,k)
             psiV(i,j+1,n,k) = 2.0*psiC-psiV(i,j,n,k)
             !  Negative flux fix up
             if (psiH(i,j,n,k)<0) then
                 psiH(i,j,n,k) = 0
             endif
             if (psiV(i,j+1,n,k)<0) then
                 psiV(i,j+1,n,k) = 0
             endif
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + w(nn)*psiC
          enddo
       enddo
       nn=nn+1
    enddo
    enddo
    !       ----------------------------------------------------------------- 
    !        mu > 0, eta > 0,   LEFT-TO-RIGHT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    psiH(1,1:TNFM_y,3*numord/4+1:numord,:) = psiH(1,1:TNFM_y,numord/2+1:3*numord/4,:) 
    psiV(1:TNFM_x,1,3*numord/4+1:numord,:) = psiV(1:TNFM_x,1,numord/4+1:numord/2,:)
    do k = 1,ng
    nn = 1
    do n = 3*numord/4+1,numord
       do j = 1,TNFM_y
             ay = 2.0*eta(nn)/dy(j)
          do i = 1,TNFM_x
             ax = 2.0*mu(nn)/dx(i)
             nume = Q_ijn(i,j,n,k) + ax*psiH(i,j,n,k) + ay*psiV(i,j,n,k)
             deno = ax + ay + SigT(fmmid2D(i,j),k)
             psiC = nume/deno
             psiH(i+1,j,n,k) = 2.0*psiC-psiH(i,j,n,k)
             psiV(i,j+1,n,k) = 2.0*psiC-psiV(i,j,n,k)
             !  Negative flux fix up
             if (psiH(i+1,j,n,k)<0) then
                 psiH(i+1,j,n,k) = 0
             endif
             if (psiV(i,j+1,n,k)<0) then
                 psiV(i,j+1,n,k) = 0
             endif
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + w(nn)*psiC
          enddo
       enddo
       nn=nn+1
    enddo
    enddo  
    phi_ij = phi
    end subroutine Sweep2D
!>BLOC10
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
!>BLOC11
    subroutine title1() 
       print*, '  ____                   _   _ _______ _____  '  
       print*, ' / __ \                 | \ | |__   __|  __ \ '
       print*, '| |  | |_ __   ___ _ __ |  \| |  | |  | |__) |'
       print*, '| |  | | "_ \ / _ \ "_ \| . ` |  | |  |  ___/ '
       print*, '| |__| | |_) |  __/ | | | |\  |  | |  | |     '
       print*, ' \____/| .__/ \___|_| |_|_| \_|  |_|  |_|     '
       print*, '       | |                                    '
       print*, '       |_|                                    '
       print*, '____________________________________________________________'
       print*, '         | The OpenNTP Neutron Transport Package            '       
       print*, 'Version  | Version Number: 1.2                              '     
       print*, 'Copyright| 2015-2020 Radiation and Nuclear system Laboratory'
       print*, '         | University Abdelmalk Essaadi, FS Tetouan, Morocco'
       print*, 'Source   | FORTRAN90 version                                ' 
       print*, 'GUI      | PyQt5                                            ' 
       print*, 'Method   | The Discrete Ordinates Method (SN)               '  
       print*, 'Dimension| Two dimensions (2D)                              '
       print*, 'Geometry | Cartesian                                        ' 
       print*, '____________________________________________________________'
       print*, ''
    end subroutine title1
!>BLOC12
    subroutine title2()
       print*,'**********************************************************'
       print*,'                         Finished                         '                             
       print*,'**********************************************************'  
    end subroutine title2
!>BLOC13
    subroutine NormalizeFlux(TNFM_x,TNFM_y,Nmat,ng,sigF,Chi,fmmid2D,xfm,yfm,phi) 
    implicit none
    ! Neutron scalar flux is normalized according to sum(V*Chi*NusigF*phi=1)
    integer(kind=4), intent(in) :: TNFM_x,TNFM_y,Nmat,ng
    real(kind=4), dimension(Nmat,ng), intent(in) :: sigF,Chi
    real(kind=4), dimension(TNFM_x), intent(in) :: xfm
    real(kind=4), dimension(TNFM_y), intent(in) :: yfm
    integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
    real(kind=4), dimension(TNFM_x,TNFM_y,ng), intent(inout) :: phi
    integer(kind=4) :: i,j
    real(kind=4) :: norme,a1
    ! Initialize local variables
    norme = 0.0
    ! Normalized source
    do j = 1,TNFM_y
        do i = 1,TNFM_x
           a1 = sum(xfm(i)*yfm(j)*sigF(fmmid2D(i,j),:)*phi(i,j,:))   
           norme = norme  + sqrt(a1*a1)
        enddo
    enddo
    phi = dot_product(xfm,yfm)*(phi/norme)
    end subroutine NormalizeFlux
!>BLOC14
    subroutine PinPower(TNFM_y,TNFM_x,Nmat,ng,nx,ny,npx,npy,npc,napc,na,nxx,nyy,xfm,yfm,&
                            core,assembly,nfmesh_xy,phi_ij,fmmid2D,SigF,PF,PD)
       ! CALCULATION OF PIN POWER 
       implicit none
       integer(kind=4), intent(in) :: TNFM_y,TNFM_x,Nmat,ng,nx,ny,npx,npy,npc,napc,na,nxx,nyy
       real(kind=4), dimension(Nmat,ng), intent(in) :: SigF
       integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
       integer(kind=4), dimension(npc,npx,npy), intent(in) ::  nfmesh_xy
       real(kind=4), dimension(TNFM_x), intent(in) :: xfm
       real(kind=4), dimension(TNFM_y), intent(in) :: yfm
       integer(kind=4), dimension(nx,ny), intent(in) :: core
       integer(kind=4), dimension(na,nxx,nyy), intent(in) :: assembly
       real(kind=4), dimension(TNFM_x,TNFM_y,ng), intent(in) :: phi_ij 
       real(kind=4), dimension(nx*nxx,ny*nyy), intent(out) :: PF,PD
       real(kind=4), dimension(nx*nxx,ny*nyy) :: PF0
       integer(kind=4) :: i,j,k,l,n1,n2,som1,som2
       real(kind=4) :: som,pm,moy
       !real(kind=4) :: VT,norme
       som1  = sum(nfmesh_xy(assembly(core(1,1),1,1),:,1)) 
       som2  = sum(nfmesh_xy(assembly(core(1,1),1,1),1,:)) 
       PF0 = 0.0; PF = 0.0; PD = 0.0

       do i = 1,nx*nxx
          n1=1
          do j = 1,ny*nyy
             do k = 1,som1
                n2 = som2*i-som2+1
                do l = 1,som2
                   PF0(i,j) = PF0(i,j) + sum(SigF(fmmid2D(n1,n2),:)*phi_ij(n1,n2,:))
                   n2 = n2+1
                enddo
                n1=n1+1
             enddo
          enddo
       enddo
       PD  = PF0
       pm  = sum(PF0)/float(napc)
       PF  = PF0/pm
    end subroutine PinPower
!>BLOC15
    subroutine AssemblyPower(nx,ny,nxx,nyy,PF,A)
       ! CALCULATION OF ASSEMBLY POWER 
       implicit none
       integer(kind=4), intent(in) :: nx,ny,nxx,nyy
       real(kind=4), dimension(nx*nxx,ny*nyy), intent(in) :: PF
       real(kind=4), dimension(nx,ny), intent(out) :: A
       integer(kind=4) :: i,j,n,m
       open (1000,file='app/Output/AssemblyPower.h')
       n=1
       do i = 1,nx
          m=1
          do j = 1,ny
             A(i,j) =  sum(PF(n:i*nxx,m:j*nyy))
             m=m+nyy
          enddo
          n=n+nxx
       enddo
       do i=1,nx
       write(1000,2000) (A(i,j),j=1,ny)
       enddo
       2000 format(10000f15.8)
    end subroutine AssemblyPower
!>BLOC16
    subroutine ScalarFluxPinC(TNFM_y,TNFM_x,Nmat,ng,nx,ny,npx,npy,npc,na,nxx,nyy,xfm,yfm,&
                            core,assembly,nfmesh_xy,phi_ij,SFPC)
    ! CALCULATION OF SCALAR FLUX IN EACH PIN CELL
       implicit none
       integer(kind=4), intent(in) :: TNFM_y,TNFM_x,Nmat,ng,nx,ny,npx,npy,npc,na,nxx,nyy
       integer(kind=4), dimension(npc,npx,npy), intent(in) ::  nfmesh_xy
       real(kind=4), dimension(TNFM_x), intent(in) :: xfm
       real(kind=4), dimension(TNFM_y), intent(in) :: yfm
       integer(kind=4), dimension(nx,ny), intent(in) :: core
       integer(kind=4), dimension(na,nxx,nyy), intent(in) :: assembly
       real(kind=4), dimension(TNFM_x,TNFM_y,ng), intent(in) :: phi_ij 
       real(kind=4), dimension(nx*nxx,ny*nyy,ng), intent(out) :: SFPC
       integer(kind=4) :: i,j,k,l,n1,n2,som1,som2

       som1  = sum(nfmesh_xy(assembly(core(1,1),1,1),:,1)) 
       som2  = sum(nfmesh_xy(assembly(core(1,1),1,1),1,:)) 
       SFPC = 0.0
       do i = 1,nx*nxx
          n1=1
          do j = 1,ny*nyy
             do k = 1,som1
                n2 = som2*i-som2+1
                do l = 1,som2
                   SFPC(i,j,:) = SFPC(i,j,:) + xfm(n1)*yfm(n2)*phi_ij(n1,n2,:)
                   n2 = n2+1
                enddo
                n1=n1+1
             enddo
          enddo
       enddo
       SFPC = SFPC/DOT_PRODUCT(xfm,yfm)
    end subroutine ScalarFluxPinC
!>BLOC17
subroutine Output(start,tm,k_eff,SigT,NusigF,SigS,Chi,mup,etap,psip,pwi,xcm,ycm,&
                  phi_ij,eps,TNFM_y,TNFM_x,ng,Nmat,order,npx,npy,N,it1,it2,npc,&
                  na,nx,ny,SFPC,PF)
        implicit none
        ! Globa Variables
        integer(kind=4), intent(in) :: ng,TNFM_y,TNFM_x,Nmat,order,npx,npy,N,it1,it2,npc,na,nx,ny
        real(kind=4), dimension(Nmat,ng), intent(in) :: SigT,NusigF,Chi
        real(kind=4), dimension(Nmat,order,ng,ng), intent(in) :: SigS
        real(kind=4), dimension(TNFM_x,TNFM_y,ng), intent(in) :: phi_ij
        real(kind=4), dimension(npc,npx), intent(in)  :: xcm
        real(kind=4), dimension(npc,npy), intent(in)  :: ycm
        real(kind=4), dimension(nx,ny,ng), intent(in) :: SFPC
        real(kind=4), dimension(nx,ny), intent(in) :: PF
        real(kind=4), dimension(N*(N+2)/8), intent(in) :: mup,etap,psip,pwi
        CHARACTER(50), intent(in) :: start,tm
        real(kind=4), intent(in) :: eps,k_eff
        ! Local Variables
        integer(kind=4) :: i,j,k
        open (100,file='app/Output/OUTPUT_SN2D.TXT')
        write (100, FMT=* ) '********************************************************************************'
        write (100, FMT=* ) 'ERSN, UNIVERSITY ABDELMALEK ESSAADI FACULTY OF SCIENCES - TETOUAN, MOROCCO'
        write (100, FMT=* ) 'CODE  DEVELOPED  BY  MOHAMED  LAHDOUR,  PHD  STUDENT'
        write (100, FMT=* ) 'OpenNTP:         SN  DISCRETE  ORDINATES  METHOD'
        write (100, FMT=* ) 'DIMENSION:       TWO DIMENSIONS (2D) '
        write (100, FMT=* ) 'GEOMETRY:        CARTESIAN'
        write (100, FMT=* ) 'VERSION NUMBER:  1.2'
        write (100, FMT=* ) 'VERSION DATE:    10  September  2020'                  
        write (100,3010) 'RAN ON:          ', start,'(H:M:S)'
        write (100, FMT=* ) '********************************************************************************'
        write (100, FMT=* ) '           ----------------------------------------------------------' 
        write (100, FMT=* ) '                     INPUT  PARAMETER - VALUES  FROM  INPUT'              
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'ENERGY GROUPS NUMBER:                    ',ng
        write (100, FMT=* ) 'PIN CELLS NUMBER:                        ',npc
        write (100, FMT=* ) 'ASSEMBLIES NUMBER:                       ',na
        write (100, FMT=* ) 'X REGIONS NUMBER PIN CELL:               ',npx
        write (100, FMT=* ) 'Y REGIONS NUMBER PIN CELL:               ',npy
        write (100, FMT=* ) 'MATERIALS NUMBER:                        ',Nmat
        do i=1,npc
        write (100,3040)    'SIZE OF EACH X REGION PIN CELL',i,':     ',xcm(i,:) 
        write (100,3040)    'SIZE OF EACH Y REGION PIN CELL',i,':     ',ycm(i,:)   
        enddo  
        write (100, FMT=* ) 'NUMBER OF DIRECTIONS ALONG EACH AXIS:    ',N
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
        write (100, FMT=* ) '--------------------------------------------------------------------------------------'
        write (100, FMT=* ) '|     N. ORDER ','          MU    ','          ETA    ','         PSI    ','       WEIGHTS       |'
        write (100, FMT=* ) '--------------------------------------------------------------------------------------'
        do i = 1,N*(N+2)/8  
          write(100,3060) '|',i,mup(i),etap(i),psip(i),pwi(i),'    |'
        enddo
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'PSEUDO  CROSS  SECTIONS  DATA: '
        write (100, FMT=* ) ''

        do i = 1,Nmat
        write (100, FMT=* ) ''
        write (100, 3070) ' MATERIAL :', i  
        write (100, FMT=* ) ''
        write (100, FMT=* ) '---------------------------------------------------------------------------&
                             &----------------------------------'
        write (100, FMT=* ) '|       GROUP ','          TOTAL ','       ABSORPTION ',&
                            '     NU*FISSION ','     SCATTERING ','     FISSION SPECTRUM       |'
        write (100, FMT=* ) '---------------------------------------------------------------------------&
                             &----------------------------------'
            do j = 1,ng
            write(100,3080) '|',j,SigT(i,j),SigT(i,j)-SigS(i,1,j,j),NusigF(i,j),SigS(i,1,j,j),Chi(i,j),'    |'
            enddo
        enddo
        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '                   REACTION RATE SOLUTION IN THE REACTOR CORE        '   
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'REACTION RATE SOLUTION PER MATERIAL PER ENERGY  GROUP: '
        write (100, FMT=* ) ''
        do i = 1,Nmat
        write (100, FMT=* ) ''
        write (100, 3070) ' MATERIAL :', i  
        write (100, FMT=* ) ''
        write (100, FMT=* ) '---------------------------------------------------------------------------&
                             &-------------'
        write (100, FMT=* ) '|       GROUP ','          TOTAL ','       ABSORPTION ',&
                            '     NU*FISSION ','     SCATTERING        |'
        write (100, FMT=* ) '---------------------------------------------------------------------------&
                             &-------------'
            do j = 1,ng
            write(100,3081) '|',j, sum(SigT(i,j)*phi_ij(:,:,j)), sum((SigT(i,j)-SigS(i,1,j,j))*phi_ij(:,:,j)),&
                               sum(NusigF(i,j)*phi_ij(:,:,j)),sum(SigS(i,1,j,j)*phi_ij(:,:,j)),'    |'
            enddo
        enddo
        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '                              NORMALIZED PIN POWER                   '   
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write(100,'(1x,2000A13)') ('-------------', i=1,nx+1)
        write(100,'(1x,A1,i14,2000i13)') '|',(i, i=1,nx)
        write(100,'(1x,2000A13)') ('-------------', i=1,nx+1)
          do j = 1,ny
          write(100,2000) '|',j,(PF(i,j), i=1,nx)
          enddo

        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '                NORMALIZED SCALAR FLUX SOLUTION IN EACH PIN CELL     '   
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
         write (100, FMT=* ) 'FLUXES  PER  PIN CELL  PER  ENERGY  GROUP:' 
        do k = 1,ng 
        write (100, FMT=* ) '' 
        write (100,3030)' P I N  C E L L','     G R O U P',k
        write (100, FMT=* ) ''
        write(100,'(1x,2000A13)') ('-------------', i=1,nx+1)
        write(100,'(1x,A1,i14,2000i13)') '|',(i, i=1,nx)
        write(100,'(1x,2000A13)') ('-------------', i=1,nx+1)
          do j = 1,ny
          write(100,2000) '|',j,(SFPC(i,j,k), i=1,nx) 
          enddo
        enddo
        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '                       NORMALIZED SCALAR  FLUX  SOLUTION             ' 
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'FLUXES  PER  MESH  PER  ENERGY  GROUP:' 
        do k = 1,ng 
        write (100, FMT=* ) '' 
        write (100,3000)' M E S H ','     G R O U P',k
        write (100, FMT=* ) ''
        write(100,'(1x,2000A13)') ('-------------', i=1,TNFM_x+1)
        write(100,'(1x,A1,i14,2000i13)') '|',(i, i=1,TNFM_x)
        write(100,'(1x,2000A13)') ('-------------', i=1,TNFM_x+1)
          do j = 1,TNFM_y
          write(100,2000) '|',j,(phi_ij(i,j,k), i=1,TNFM_x)
          enddo
        enddo
        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '             OUTPUT  PARAMETER - SOLUTION  TO  TRANSPORT  EQUATION   ' 
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100,3090)    'K-EFF                    =',k_eff
        write (100,3020)    'TOTAL OUTER ITERATIONS   =',it1
        write (100,3020)    'TOTAL INNER ITERATIONS   =',it2
        write (100,4000)    'TOTAL EXECUTION TIME     =',tm,'(H:M:S)'
        write (100, FMT=* ) ''
        write (100, FMT=* ) '********************************************************************************'
        2000 format(1x,1p,A1,i4,1x,1000000e13.5) 
        3000 format(1x,A8,2x,300(A14,i2))  
        3010 format(1x,A17,A22,A10)
        3020 format(1x,A26,4x,i10)
        3030 format(1x,A18,2x,300(A14,i2))
        3040 format(1x,A30,i2,2x,A6,9x,200F10.5)
        3050 format(1x,1p,A41,4x,e8.1)
        3060 format(1x,1p,A1,i11,5x,e16.5,e16.5,e16.5,e16.5,A5)
        3070 format(1x,A18,i4)
        3080 format(1x,1p,A1,i11,5x,e16.5,e16.5,e16.5,e16.5,e16.5,A12)
        3081 format(1x,1p,A1,i11,5x,e16.5,e16.5,e16.5,e16.5,A7)
        3090 format(1x,A26,6x,f8.6)
        4000 format(1x,A26,4x,A10,A10)
        close(100)
end subroutine Output      
