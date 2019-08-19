!BL1
    subroutine quadrature_set(N,mup, etap, psip,pwi)
       !Calculate the quadrature sets in an octant and return
       !complete angle and weight set over all
       implicit none
       integer(kind=8), intent(in) :: N
       real(kind=8), dimension(N/2) :: mmu
       real(kind=8), dimension(N*(N+2)/8), intent(out) :: mup, etap, psip,pwi
       real(kind=8), parameter :: pi = 3.141592654
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
    subroutine fmm_id(fmmid_xy,nfmesh_xy,RegMat,nreg_x,nreg_y,totnfm_y,totnfm_x)
       implicit none
       integer(kind=4), intent(in) :: nreg_x,nreg_y,totnfm_y,totnfm_x
       integer(kind=4), dimension(nreg_y,nreg_x), intent(in) :: RegMat
       integer(kind=4), dimension(nreg_y,nreg_x), intent(in) :: nfmesh_xy
       integer(kind=4), dimension(totNFM_y,totnfm_x), intent(out) :: fmmid_xy
       integer(kind=4) :: i,j,m,n,nn,mm,k0
       !-- fine mesh material id x direction     
       print*,'hi man'
       k0 = 1
       do i=1,nreg_x
       mm = 1
          do j=1,nreg_y             
             do m = 1,nfmesh_xy(j,i) 
                nn=k0 
                do n = 1,nfmesh_xy(j,i)
                   fmmid_xy(mm,nn) = RegMat(j,i)
                   nn=nn+1
                enddo
                mm=mm+1
             enddo
          enddo
          k0=nn
       enddo 
    end subroutine fmm_id
!BL3 
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
!BL4
    subroutine Plot_plm(l,m,A)
        implicit none
        ! Variable Globale
        integer(kind=4),intent(in) :: l,m
        real(kind=8), dimension(l+1,m+1,101),intent(out) :: A
        ! Variable Locale
        real(kind=8), dimension(101) :: mu
        real(kind=8) :: x,plm 
        integer(kind=4) :: i,j,k
        Open (10,file='app/Output/plot_plm.h')
        x = 0.
        mu = 0.
        do i = 1,101
        mu(i) = x 
               x = x + 0.00999995
        enddo
        do i=1,l+1         ! itération sur m
           do j=i,m+1      ! itération sur l
              do k=1,101
                 A(i,j,k)=plm(j-1,i-1,mu(k))
              enddo
           enddo
        enddo
        x = 0.
        do k=1,101
        write(10,'(1000f10.6)') x,((A(j,i,k),i=j,m+1), j=1,l+1)
        x = x + 0.00999995
        enddo
    end subroutine Plot_plm  
!BL5    
    subroutine Delta_f(nreg_x,nreg_y,totNFM_x,totNFM_y,nfmesh_xy,cell_x,cell_y,Delta_x,Delta_y)
       implicit none
       integer(kind=4), intent(in) :: nreg_x,nreg_y,totNFM_x,totNFM_y
       integer(kind=4), dimension(nreg_y,nreg_x), intent(in) :: nfmesh_xy
       real(kind=8), dimension(nreg_x), intent(in) :: cell_x
       real(kind=8), dimension(nreg_y), intent(in) :: cell_y
       real(kind=8), dimension(totNFM_x,nreg_y), intent(out) :: Delta_x
       real(kind=8), dimension(nreg_x,totNFM_y), intent(out) :: Delta_y
       integer(kind=4) :: i,j,k,m,n
       !-- Size each mesh in the x direction 
       m = 1
       do j = 1,nreg_y
          n = 1
          do i = 1,nreg_x

              do k = 1,nfmesh_xy(j,i)
                 Delta_x(n,m) = cell_x(i)/nfmesh_xy(j,i)
                 n = n+1
              enddo
          enddo
          m = m+1
       enddo
       !-- Size each mesh in the y direction
       m = 1
       do i = 1,nreg_x
          n = 1
          do j = 1,nreg_y
              do k = 1,nfmesh_xy(j,i)
                 Delta_y(m,n) = cell_y(i)/nfmesh_xy(j,i)
                 n = n+1
              enddo
          enddo
          m = m+1
       enddo
    end subroutine Delta_f
!BL6
    subroutine Matrix_AxAyB(ng,Nmat,totNFM_x,totNFM_y,ngauss,nreg_y,nreg_x,&
                            mup,etap,fmmid_xy,SigT,Delta_x,Delta_y,B,Ax,Ay)
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,totNFM_x,totNFM_y,ngauss,nreg_y,nreg_x
       real(kind=8), dimension(totNFM_x,nreg_y), intent(in) :: Delta_x
       real(kind=8), dimension(nreg_x,totNFM_y), intent(in) :: Delta_y
       integer(kind=4), dimension(totNFM_y,totnfm_x), intent(in) :: fmmid_xy
       real(kind=8), dimension(Nmat,ng), intent(in) :: SigT
       real(kind=8), dimension(ngauss*(ngauss+2)/8), intent(in) :: mup, etap
       real(kind=8), dimension(totNFM_x,totnfm_y,ngauss*(ngauss+2)/8), intent(out) :: B,Ax,Ay
       integer(kind=4) :: i,j,n
       do n = 1,ngauss*(ngauss+2)/8
          do j = 1,totNFM_y
             do i = 1,totNFM_x
                Ax(i,j,n) = 2.d0*mup(n)/Delta_x(i,fmmid_xy(1,i))
                Ay(i,j,n) = 2.d0*etap(n)/Delta_y(fmmid_xy(j,1),j)
                B(i,j,n)  = SigT(fmmid_xy(j,i),1) + Ax(i,j,n) + Ay(i,j,n)
             enddo
          enddo
       enddo
    end subroutine Matrix_AxAyB
!BL7
    subroutine Flux_guess(phi_lmij,totNFM_x,totNFM_y,ngauss,mup,pwi,order)
        implicit none
        integer(kind=4), intent(in) :: totNFM_x,totNFM_y,ngauss,order
        real(kind=8), dimension(ngauss*(ngauss+2)/8), intent(in) :: mup,pwi
        real(kind=8), dimension(totNFM_x,totnfm_y,order,order), intent(out) :: phi_lmij
        real(kind=8), dimension(totNFM_x,totnfm_y,ngauss*(ngauss+2)/2) :: phi_ijn
        real(kind=8) :: plm
        integer(kind=4) :: i,j,ii,n,nn,l,m
        phi_ijn = 1      
        do l=0,order-1     ! itération sur l
        do m=0,l           ! itération sur m   
        do j = 1,totNFM_y
           do i = 1,totNFM_x
              nn=1
              do ii =1,4
                 do n = 1,ngauss*(ngauss+2)/8
                    phi_lmij(i,j,l+1,m+1) = phi_lmij(i,j,l+1,m+1) + &
                         0.25D0*pwi(n)*phi_ijn(i,j,nn)*plm(l,m,mup(n))
                    nn=nn+1
                 enddo 
              enddo
           enddo
        enddo
        enddo
        enddo
    end subroutine Flux_guess
!BL8
    subroutine Phi_lmij_f(phi_lmij,phi_ijn,totNFM_x,totNFM_y,ngauss,mup,pwi,order)
        implicit none
        integer(kind=4), intent(in) :: totNFM_x,totNFM_y,ngauss,order
        real(kind=8), dimension(ngauss*(ngauss+2)/8), intent(in) :: mup,pwi
        real(kind=8), dimension(totNFM_x,totnfm_y,ngauss*(ngauss+2)/2), intent(in) :: phi_ijn
        real(kind=8), dimension(totNFM_x,totnfm_y,order,order), intent(out) :: phi_lmij
        real(kind=8) :: plm
        integer(kind=4) :: i,j,ii,n,nn,l,m    
        phi_lmij = 0.0D0
        do l=0,order-1     ! itération sur l
        do m=0,l           ! itération sur m   
        do j = 1,totNFM_y
           do i = 1,totNFM_x
              nn=1
              do ii =1,4
                 do n = 1,ngauss*(ngauss+2)/8
                    phi_lmij(i,j,l+1,m+1) = phi_lmij(i,j,l+1,m+1) + &
                         0.25D0*pwi(n)*phi_ijn(i,j,nn)*plm(l,m,mup(n))
                    nn=nn+1
                 enddo 
              enddo
           enddo
        enddo
        enddo
        enddo
    end subroutine Phi_lmij_f
!BL9
    subroutine Scattering_Source(totNFM_x,totNFM_y,Nmat,order,ng,ngauss,&
                                 fmmid_xy,mup,SigS,phi_lmij,Qs_nij)
        implicit none
        integer(kind=4), intent(in) :: totNFM_x,totNFM_y,Nmat,order,ng,ngauss
        real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
        real(kind=8), dimension(ngauss*(ngauss+2)/8), intent(in) :: mup
        integer(kind=4), dimension(totNFM_y,totnfm_x), intent(in) :: fmmid_xy
        real(kind=8), dimension(totNFM_x,totnfm_y,order,order), intent(in) :: phi_lmij
        real(kind=8), dimension(totNFM_x,totnfm_y,ngauss*(ngauss+2)/2), intent(out) :: Qs_nij
        real(kind=8), dimension(totNFM_x,totnfm_y,order,order) :: Qs_lmij
        real(kind=8), parameter :: PI = 3.141592654
        real(kind=8) :: plm,som
        integer(kind=4) :: i,j,l,m,nn,n,ii
        Qs_nij = 0.0D0
        do l = 0,order-1
           do m = 0,l
              do j = 1,totNFM_y
                 do i = 1,totNFM_x
                    Qs_lmij(i,j,l+1,m+1) = SigS(fmmid_xy(j,i),l+1,1,1)*phi_lmij(i,j,l+1,m+1)
                 enddo
              enddo
           enddo
        enddo

        do j = 1,totNFM_y
           do i = 1,totNFM_x
              som = 0.0D0
              do l = 0,order-1
              nn = 1
              do ii =1,4
                 do n = 1,ngauss*(ngauss+2)/8
                    do m = 0,l
                       som = som + Qs_lmij(i,j,l+1,m+1)*plm(l,m,mup(n))
                    enddo
                    Qs_nij(i,j,nn) = Qs_nij(i,j,nn) + (0.25D0/PI)*DBLE(2*l+1)*som
                    nn=nn+1
                 enddo
              enddo
              enddo
           enddo
        enddo
    end subroutine Scattering_Source
!BL10
    subroutine Fission_Source(totNFM_x,totNFM_y,Nmat,order,ng,ngauss,fmmid_xy,&
                                          mup,NusigF,Chi,phi_lmij,k_eff,Qf_nij)
        implicit none
        integer(kind=4), intent(in) :: totNFM_x,totNFM_y,Nmat,order,ng,ngauss
        real(kind=8), intent(in) :: k_eff 
        real(kind=8), dimension(Nmat,ng), intent(in) :: NusigF,Chi
        real(kind=8), dimension(ngauss*(ngauss+2)/8), intent(in) :: mup
        integer(kind=4), dimension(totNFM_y,totnfm_x), intent(in) :: fmmid_xy
        real(kind=8), dimension(totNFM_x,totnfm_y,order,order), intent(in) :: phi_lmij
        real(kind=8), dimension(totNFM_x,totnfm_y,ngauss*(ngauss+2)/2), intent(out) :: Qf_nij
        real(kind=8), dimension(totNFM_x,totnfm_y,order,order) :: Qf_lmij,phi_11ij
        real(kind=8), parameter :: PI = 3.141592654
        real(kind=8) :: plm,som
        integer(kind=4) :: i,j,l,m,nn,n,ii
        Qf_lmij = 0.0D0
        phi_11ij = 0
        phi_11ij(:,:,1,1) = phi_lmij(:,:,1,1)

        do l = 0,order-1
           do m = 0,l
              do j = 1,totNFM_y
                 do i = 1,totNFM_x
                    Qf_lmij(i,j,l+1,m+1) = Chi(fmmid_xy(j,i),1)*&
                           NusigF(fmmid_xy(j,i),1)*phi_11ij(i,j,l+1,m+1)/k_eff
                 enddo
              enddo
           enddo
        enddo


        do j = 1,totNFM_y
           do i = 1,totNFM_x
              som = 0.0D0
              do l = 0,order-1
              nn = 1
              do ii =1,4
                 do n = 1,ngauss*(ngauss+2)/8
                    do m = 0,l
                       som = som + Qf_lmij(i,j,l+1,m+1)*plm(l,m,mup(n))
                    enddo
                    Qf_nij(i,j,nn) = Qf_nij(i,j,nn) + (0.25D0/PI)*DBLE(2*l+1)*som
                    nn=nn+1
                 enddo
              enddo
              enddo
           enddo
        enddo
    end subroutine Fission_Source
!BL11
    subroutine Sweeping(BC1,BC2,BC3,BC4,totNFM_x,totnfm_y,ngauss,Q_nij,B,Ax,Ay,phi_ijn)
    implicit none
    CHARACTER(50), intent(in) :: BC1,BC2,BC3,BC4
    integer(kind=4), intent(in) :: totNFM_x,totnfm_y,ngauss
    real(kind=8), dimension(totNFM_x,totnfm_y,ngauss*(ngauss+2)/2), intent(in) :: Q_nij 
    real(kind=8), dimension(totNFM_x,totnfm_y,ngauss*(ngauss+2)/8), intent(in) :: B,Ax,Ay
    real(kind=8), dimension(totNFM_x,totnfm_y,ngauss*(ngauss+2)/2), intent(out) :: phi_ijn
    real(kind=8), dimension(totNFM_x+1,totnfm_y,ngauss*(ngauss+2)/8) :: currx,currax,currbx
    real(kind=8), dimension(totNFM_x,totnfm_y+1,ngauss*(ngauss+2)/8) :: curry,curray,currby
    integer(kind=4) :: i,j,n,nn

    if (BC1 =='Reflective Right') then
    !========================\\reflective in the right side//==========================
    ! Left to right; bottom to top
    do n = 1,ngauss*(ngauss+2)/8
       do j = 1,totNFM_y
          currax(1,j,n) = 0.0D0 ! pour mu>0
          currbx(1,j,n) = 1.0D0
          do i = 1,totNFM_x         
             curray(i,1,n)   = 0.0D0
             currby(i,1,n)   = 1.0D0
             phi_ijn(i,j,n) = (Ax(i,j,n)*currax(i,j,n)+Ay(i,j,n)*curray(i,j,n)+Q_nij(i,j,n))/B(i,j,n)
             phi_ijn(i,j,n) = (Ax(i,j,n)*currbx(i,j,n)+Ay(i,j,n)*currby(i,j,n)+Q_nij(i,j,n))/B(i,j,n)
             currax(i+1,j,n) = 2.0D0*phi_ijn(i,j,n)-currax(i,j,n)
             currbx(i+1,j,n) = 2.0D0*phi_ijn(i,j,n)-currbx(i,j,n)
             currx(i+1,j,n)= currax(i,j,n)/ (1.0D0 + (currax(i,j,n)-currbx(i,j,n)))
          enddo
          do i = 1,totNFM_x
             curray(i,j+1,n) = 2.0D0*phi_ijn(i,j,n)-curray(i,j,n)
             currby(i,j+1,n) = 2.0D0*phi_ijn(i,j,n)-currby(i,j,n)
             curry(i,j+1,n)= curray(i,j,n)/ (1.0D0 + (curray(i,j,n)-currby(i,j,n)))
          enddo
       enddo
    enddo
    do n = 1,ngauss*(ngauss+2)/8
       do j = 1,totNFM_y
          do i = 1,totNFM_x         
             phi_ijn(i,j,n) = (Ax(i,j,n)*currx(i,j,n)+Ay(i,j,n)*curry(i,j,n)+Q_nij(i,j,n))/B(i,j,n)
             currx(i+1,j,n) = 2.0D0*phi_ijn(i,j,n)-currx(i,j,n)
          enddo
          do i = 1,totNFM_x
             curry(i,j+1,n) = 2.0D0*phi_ijn(i,j,n)-curry(i,j,n)
          enddo
       enddo
    enddo
    else
    !========================\\vacum in the right side//==========================
    ! Left to right; bottom to top
    do n = 1,ngauss*(ngauss+2)/8
       do j = 1,totNFM_y
          currx(1,j,n) = 0.0D0 ! pour mu>0
          do i = 1,totNFM_x         
             curry(i,1,n)   = 0.0D0
             phi_ijn(i,j,n) = (Ax(i,j,n)*currx(i,j,n)+Ay(i,j,n)*curry(i,j,n)+Q_nij(i,j,n))/B(i,j,n)
             currx(i+1,j,n) = 2.0D0*phi_ijn(i,j,n)-currx(i,j,n)
          enddo
          do i = 1,totNFM_x
             curry(i,j+1,n) = 2.0D0*phi_ijn(i,j,n)-curry(i,j,n)
          enddo
       enddo
    enddo
    endif


    if (BC2 == 'Reflective Left') then
    !========================\\reflective in the left side//==========================
    ! Right to left; bottom to top
    nn = n
    do n = 1,ngauss*(ngauss+2)/8
       do j = 1,totNFM_y
          currax(totNFM_x+1,j,n) = 0.0D0 ! pour mu<0
          currbx(totNFM_x+1,j,n) = 1.0D0
          do i = totNFM_x,1,-1        
             curray(i,1,n) = 0.0D0
             currby(i,1,n) = 1.0D0
             phi_ijn(i,j,nn) = (Ax(i,j,n)*currax(i+1,j,n)+Ay(i,j,n)*curray(i,j,n)+Q_nij(i,j,n))/B(i,j,n)
             phi_ijn(i,j,nn) = (Ax(i,j,n)*currbx(i+1,j,n)+Ay(i,j,n)*currby(i,j,n)+Q_nij(i,j,n))/B(i,j,n)
             currax(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-currax(i+1,j,n)
             currbx(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-currbx(i+1,j,n)
             currx(i,j,n)= currax(i+1,j,n)/ (1.0D0 + (currax(i+1,j,n)-currbx(i+1,j,n)))
          enddo
          do i = totNFM_x,1,-1
             curray(i,j+1,n) = 2.0D0*phi_ijn(i,j,nn)-curray(i,j,n)
             currby(i,j+1,n) = 2.0D0*phi_ijn(i,j,nn)-currby(i,j,n)
             curry(i,j+1,n)= curray(i,j,n)/ (1.0D0 + (curray(i,j,n)-currby(i,j,n)))
          enddo
       enddo
    nn=nn+1
    enddo
    nn = n
    do n = 1,ngauss*(ngauss+2)/8
       do j = 1,totNFM_y
          do i = totNFM_x,1,-1        
             phi_ijn(i,j,nn) = (Ax(i,j,n)*currx(i+1,j,n)+Ay(i,j,n)*curry(i,j,n)+Q_nij(i,j,n))/B(i,j,n)
             currx(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-currx(i+1,j,n)
          enddo
          do i = totNFM_x,1,-1
             curry(i,j+1,n) = 2.0D0*phi_ijn(i,j,nn)-curry(i,j,n)
          enddo
       enddo
    nn=nn+1
    enddo
    else
    !========================\\vacuum in the left side//==========================
    ! Right to left; bottom to top
    nn = n
    do n = 1,ngauss*(ngauss+2)/8
       do j = 1,totNFM_y
          currx(totNFM_x+1,j,n) = 0.0D0 ! pour mu<0
          do i = totNFM_x,1,-1        
             curry(i,1,n) = 0.0D0
             phi_ijn(i,j,nn) = (Ax(i,j,n)*currx(i+1,j,n)+Ay(i,j,n)*curry(i,j,n)+Q_nij(i,j,n))/B(i,j,n)
             currx(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-currx(i+1,j,n)
          enddo
          do i = totNFM_x,1,-1
             curry(i,j+1,n) = 2.0D0*phi_ijn(i,j,nn)-curry(i,j,n)
          enddo
       enddo
    nn=nn+1
    enddo
    endif


    if (BC3 == 'Reflective Bottom') then
    !========================\\reflective in the bottom side//==========================
    ! Left to right; top to bottom
    do n = 1,ngauss*(ngauss+2)/8
       do j = totNFM_y,1,-1
          currax(1,j,n) = 0.0D0 ! pour mu>0
          currbx(1,j,n) = 1.0D0 
          do i = 1,totNFM_x         
             curray(i,totNFM_y+1,n) = 0.0D0
             currby(i,totNFM_y+1,n) = 1.0D0
             phi_ijn(i,j,nn)  = (Ax(i,j,n)*currax(i,j,n)+Ay(i,j,n)*curray(i,j+1,n)+Q_nij(i,j,n))/B(i,j,n)
             phi_ijn(i,j,nn)  = (Ax(i,j,n)*currbx(i,j,n)+Ay(i,j,n)*currby(i,j+1,n)+Q_nij(i,j,n))/B(i,j,n)
             currax(i+1,j,n) = 2.0D0*phi_ijn(i,j,nn)-currax(i,j,n)
             currbx(i+1,j,n) = 2.0D0*phi_ijn(i,j,nn)-currbx(i,j,n)
             currx(i+1,j,n)= currax(i,j,n)/ (1.0D0 + (currax(i,j,n)-currbx(i,j,n)))
          enddo
          do i = 1,totNFM_x
             curray(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-curray(i,j+1,n)
             currby(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-currby(i,j+1,n)
             curry(i,j,n)= curray(i,j+1,n)/ (1.0D0 + (curray(i,j+1,n)-currby(i,j+1,n)))
          enddo
       enddo
    nn=nn+1
    enddo
    nn=nn-ngauss*(ngauss+2)/8
    do n = 1,ngauss*(ngauss+2)/8
       do j = totNFM_y,1,-1
          do i = 1,totNFM_x         
             phi_ijn(i,j,nn)  = (Ax(i,j,n)*currx(i,j,n)+Ay(i,j,n)*curry(i,j+1,n)+Q_nij(i,j,n))/B(i,j,n)
             currx(i+1,j,n) = 2.0D0*phi_ijn(i,j,nn)-currx(i,j,n)
          enddo
          do i = 1,totNFM_x
             curry(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-curry(i,j+1,n)
          enddo
       enddo
    nn=nn+1
    enddo
    else
    !========================\\vacuum in the bottom side//==========================
    ! Left to right; top to bottom
    do n = 1,ngauss*(ngauss+2)/8
       do j = totNFM_y,1,-1
          currx(1,j,n) = 0.0D0 ! pour mu>0
          do i = 1,totNFM_x         
             curry(i,totNFM_y+1,n) = 0.0D0
             phi_ijn(i,j,nn)  = (Ax(i,j,n)*currx(i,j,n)+Ay(i,j,n)*curry(i,j+1,n)+Q_nij(i,j,n))/B(i,j,n)
             currx(i+1,j,n) = 2.0D0*phi_ijn(i,j,nn)-currx(i,j,n)
          enddo
          do i = 1,totNFM_x
             curry(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-curry(i,j+1,n)
          enddo
       enddo
    nn=nn+1
    enddo
    endif


    if (BC4 == 'Reflective Top') then
    !========================\\reflective in the top side//==========================
    ! Right to left; top to bottom
    do n = 1,ngauss*(ngauss+2)/8
       do j = totNFM_y,1,-1
          currax(totNFM_x+1,j,n) = 0.0D0 ! pour mu<0
          currbx(totNFM_x+1,j,n) = 1.0D0 ! pour mu<0
          do i = totNFM_x,1,-1         
             curray(i,totNFM_y+1,n) = 0.0D0
             currby(i,totNFM_y+1,n) = 1.0D0
             phi_ijn(i,j,nn)  = (Ax(i,j,n)*currax(i+1,j,n)+Ay(i,j,n)*curray(i,j+1,n)+Q_nij(i,j,n))/B(i,j,n)
             phi_ijn(i,j,nn)  = (Ax(i,j,n)*currbx(i+1,j,n)+Ay(i,j,n)*currby(i,j+1,n)+Q_nij(i,j,n))/B(i,j,n)
             currax(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-currax(i+1,j,n)
             currbx(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-currbx(i+1,j,n)
             currx(i,j,n)= currax(i+1,j,n)/ (1.0D0 + (currax(i+1,j,n)-currbx(i+1,j,n)))
          enddo
          do i = totNFM_x,1,-1 
             curray(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-curray(i,j+1,n)
             currby(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-currby(i,j+1,n)
             curry(i,j,n)= curray(i,j+1,n)/ (1.0D0 + (curray(i,j+1,n)-currby(i,j+1,n)))
          enddo
       enddo
    nn=nn+1
    enddo
    nn = nn-ngauss*(ngauss+2)/8
    do n = 1,ngauss*(ngauss+2)/8
       do j = totNFM_y,1,-1
          do i = totNFM_x,1,-1         
             phi_ijn(i,j,nn)  = (Ax(i,j,n)*currx(i+1,j,n)+Ay(i,j,n)*curry(i,j+1,n)+Q_nij(i,j,n))/B(i,j,n)
             currx(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-currx(i+1,j,n)
          enddo
          do i = totNFM_x,1,-1 
             curry(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-curry(i,j+1,n)
          enddo
       enddo
    nn=nn+1
    enddo
    else
    !========================\\vacuum in the top side//==========================
    ! Right to left; top to bottom
    do n = 1,ngauss*(ngauss+2)/8
       do j = totNFM_y,1,-1
          currx(totNFM_x+1,j,n) = 0.0D0 ! pour mu<0
          do i = totNFM_x,1,-1         
             curry(i,totNFM_y+1,n) = 0.0D0
             phi_ijn(i,j,nn)  = (Ax(i,j,n)*currx(i+1,j,n)+Ay(i,j,n)*curry(i,j+1,n)+Q_nij(i,j,n))/B(i,j,n)
             currx(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-currx(i+1,j,n)
          enddo
          do i = totNFM_x,1,-1 
             curry(i,j,n) = 2.0D0*phi_ijn(i,j,nn)-curry(i,j+1,n)
          enddo
       enddo
    nn=nn+1
    enddo
    endif
    end subroutine Sweeping
!BL12
    subroutine ScalarFlux(totNFM_x,totnfm_y,ngauss,phi_ijn,pwi,phi_ij)
    implicit none
    integer(kind=4), intent(in) :: totNFM_x,totnfm_y,ngauss
    real(kind=8), dimension(totNFM_x,totnfm_y,ngauss*(ngauss+2)/2), intent(in) :: phi_ijn
    real(kind=8), dimension(ngauss*(ngauss+2)/8), intent(in) :: pwi
    real(kind=8), dimension(totNFM_x,totnfm_y), intent(out) :: phi_ij
    real(kind=8), dimension(totNFM_x,totnfm_y)  :: phi 
    integer(kind=4) :: i,j,n,nn,ii
    phi = 0.0
        do j = 1,totNFM_y
           do i = 1,totNFM_x
              nn=1
              do ii =1,4
                 do n = 1,ngauss*(ngauss+2)/8
                    phi(i,j) = phi(i,j) + 0.25*pwi(n)*phi_ijn(i,j,nn)
                    nn=nn+1
                 enddo 
              enddo
           enddo
        enddo
    phi_ij = phi
    end subroutine ScalarFlux
!BL13
    subroutine plot_flux(totNFM_x,totnfm_y,nreg_y,nreg_x,Delta_x,Delta_y,fmmid_xy,phi_ij)
       implicit none
       integer(kind=4), intent(in) :: totNFM_x,totnfm_y,nreg_y,nreg_x
       real(kind=8), dimension(totNFM_x,nreg_y), intent(in) :: Delta_x
       real(kind=8), dimension(nreg_x,totNFM_y), intent(in) :: Delta_y
       integer(kind=4), dimension(totNFM_y,totnfm_x), intent(in) :: fmmid_xy
       real(kind=8), dimension(totNFM_x,totnfm_y), intent(in) :: phi_ij
       real(kind=8) :: somx,somy,r=0
       integer(kind=4) :: i,j
       open (14,file='app/Output/flux_sn2D.h') 
       somx =  Delta_x(1,1)
       somy =  Delta_y(1,1)

       write(14,'(f13.8,1000f13.8)') r,(sum(Delta_y(fmmid_xy(1:i,1),i)), i=1,totnfm_y) 
       do i=1,totNFM_x
          write(14,'(1000f13.8)') somx, ( phi_ij(i,j), j=1,totnfm_y)  
          somx = somx + Delta_x(i,fmmid_xy(1,i))
       enddo
       close(14)
    end subroutine plot_flux
!BL14
    subroutine Eigenvalues(totNFM_x,totNFM_y,Nmat,order,ng,N,nreg_y,nreg_x,&
                           fmmid_xy,mup,etap,pwi,Delta_x,Delta_y,SigT,NusigF,&
                           BC1,BC2,BC3,BC4,Chi,SigS,phi_ij)
    implicit none
    ! Variables Globales
    CHARACTER(50), intent(in) :: BC1,BC2,BC3,BC4
    integer(kind=4), intent(in) :: totNFM_x,totNFM_y,Nmat,order,ng,N,nreg_y,nreg_x
    integer(kind=4), dimension(totNFM_y,totnfm_x), intent(in) :: fmmid_xy
    real(kind=8), dimension(N*(N+2)/8), intent(in) :: mup,etap,pwi
    real(kind=8), dimension(totNFM_x,nreg_y), intent(in) :: Delta_x
    real(kind=8), dimension(nreg_x,totNFM_y), intent(in) :: Delta_y
    real(kind=8), dimension(Nmat,ng), intent(in) :: SigT,NusigF,Chi
    real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
    real(kind=8), dimension(totNFM_x,totnfm_y), intent(out) :: phi_ij
    ! Variables Locales
    real(kind=8), dimension(totNFM_x,totnfm_y,N*(N+2)/2) :: phi_ijn,Qf_nij,Qs_nij,Q_nij,phi_ijn0
    real(kind=8), dimension(totNFM_x,totnfm_y) :: phi_ij0,F
    real(kind=8), dimension(totNFM_x,totnfm_y,order,order) :: phi_lmij
    real(kind=8), dimension(totNFM_x,totnfm_y,N*(N+2)/8) :: B,Ax,Ay
    real(kind=8), dimension(totNFM_x*totnfm_y) :: Vect
    real(kind=8) :: err_k_eff, err_phi, eps=1.0E-8,k_eff,k_eff0,err_flux,eval1,eval2,moy
    integer(kind=4) :: inter,exter,i,j,k

    call Matrix_AxAyB(ng,Nmat,totNFM_x,totNFM_y,N,nreg_y,nreg_x,&
                            mup,etap,fmmid_xy,SigT,Delta_x,Delta_y,B,Ax,Ay)
    call Flux_guess(phi_lmij,totNFM_x,totNFM_y,N,mup,pwi,order)
    phi_ij0   = phi_lmij(:,:,1,1)
    inter     = 0
    err_k_eff = 1.0D0
    err_phi   = 1.0D0
    err_flux  = 1.0D0
    exter     = 0 
    k_eff     = 1.0D0 ! Given K_eff 
    phi_ijn0    = 1.0D0
    ! Staring extern Iteration
    do while ( err_k_eff >= eps .and. err_phi >= eps )
       call Fission_Source(totNFM_x,totNFM_y,Nmat,order,ng,N,fmmid_xy,mup,&
                           NusigF,Chi,phi_lmij,k_eff,Qf_nij)
       k_eff0   = k_eff
       err_flux = 1.0D0
       ! Starting intern Iteration
       do while (  err_flux >= eps )
          call Scattering_Source(totNFM_x,totNFM_y,Nmat,order,ng,N,fmmid_xy,&
                                 mup,SigS,phi_lmij,Qs_nij)
          Q_nij = Qf_nij + Qs_nij
          call Sweeping(BC1,BC2,BC3,BC4,totNFM_x,totnfm_y,N,Q_nij,B,Ax,Ay,phi_ijn)
          call Phi_lmij_f(phi_lmij,phi_ijn,totNFM_x,totNFM_y,N,mup,pwi,order)
          ! Condition sur le flux angulaire
          err_flux = maxval(abs(phi_ijn-phi_ijn0)/phi_ijn0)
          phi_ijn0 = phi_ijn
          inter = inter + 1
          ! ending intern Iteration 
       enddo
       phi_ij = phi_lmij(:,:,1,1)
       
       ! Exprésiion de Normalisation

       ! Condition sur le flux scalaire
       err_phi =  maxval(abs(phi_ij-phi_ij0)/phi_ij0)
       exter  = exter + 1
       eval2=0
       do j = 1,totNFM_y
           do i = 1,totNFM_x
              eval2 = eval2  + Chi(fmmid_xy(j,i),1)*NusigF(fmmid_xy(j,i),1)*phi_ij(j,i)
              eval1 = eval1  + Chi(fmmid_xy(j,i),1)*NusigF(fmmid_xy(j,i),1)*phi_ij0(j,i)
           enddo
       enddo
       k_eff = k_eff*(eval2/eval1) 
       eval1=eval2
       err_k_eff =  abs(k_eff-k_eff0)/k_eff
       write(*,2000)exter,k_eff,err_k_eff
       phi_ij0 = phi_ij
       !call Phi_lmij_f(phi_lmij,phi_ijn,totNFM_x,totNFM_y,N,mup,pwi,order)
       ! ending extern Iteration
       enddo
       call ScalarFlux(totNFM_x,totnfm_y,N,phi_ijn,pwi,phi_ij)
       k = 1
       do j = 1,totNFM_y
          do i = 1,totNFM_x
             Vect(k) = phi_ij(j,i)
             k=k+1
          enddo
       enddo
       moy = sum(Vect)/(totnfm_y*totnfm_x)
       phi_ij = phi_ij/moy
       phi_ij=phi_ij/maxval(phi_ij)
       2000 format("Iteration",i4,":",5x,"===>",5x,"keff =",F9.6,5x,"===>",5x,"res =",e10.3)
    end subroutine Eigenvalues









