     program test
     implicit none
     integer(kind=4) :: TNFM_x,TNFM_y,Nmat,order,ng,NORD
        real(kind=8), dimension(Nmat,order,ng,ng) :: SigS
        real(kind=8), dimension(NORD*(NORD+2)/8) :: mu,w
        integer(kind=4), dimension(TNFM_x,TNFM_y) :: fmmid2D
        real(kind=8), dimension(TNFM_x,TNFM_y,ng) :: phi_ij
        real(kind=8), dimension(TNFM_x,TNFM_y,ng):: Qs_ij
     call Scattering_Source(TNFM_x,TNFM_y,Nmat,order,ng,NORD,&
                                 fmmid2D,mu,w,SigS,phi_ij,Qs_ij)
     end program
     subroutine Scattering_Source(TNFM_x,TNFM_y,Nmat,order,ng,NORD,&
                                 fmmid2D,mu,w,SigS,phi_ij,Qs_ij)
        implicit none
        integer(kind=4), intent(in) :: TNFM_x,TNFM_y,Nmat,order,ng,NORD
        real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
        real(kind=8), dimension(NORD*(NORD+2)/8), intent(in) :: mu,w
        integer(kind=4), dimension(TNFM_x,TNFM_y), intent(in) :: fmmid2D
        real(kind=8), dimension(TNFM_x,TNFM_y,ng), intent(in) :: phi_ij
        real(kind=8), dimension(TNFM_x,TNFM_y,ng), intent(out) :: Qs_ij
        ! Variables locales
        real(kind=8), dimension(TNFM_x,TNFM_y,order,order,ng) :: phi_lmij
        real(kind=8), dimension(TNFM_x,TNFM_y,order,order,ng) :: Qs_lmij
        real(kind=8), dimension(TNFM_x,TNFM_y,NORD*(NORD+2)/2,ng) :: Qs_nij
        real(kind=8), dimension(TNFM_x,TNFM_y,ng) :: Qs
        real(kind=8), parameter :: PI = 3.141592653589793
        real(kind=8) :: plm,som
        integer(kind=4) :: i,j,k,l,m,n,n1,i1,numord
        numord=NORD*(NORD+2)/8
        phi_lmij(:,:,1,1,:) = phi_ij
        !>--------------- 
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
        !>--------------- 
        
        do k = 1,ng
           do j = 1,TNFM_y
           do i = 1,TNFM_x
              n1 = 1
              do i1 = 1,4
               do n  = 1,numord
               do l = 0,order-1
                     som = 0.0D0
                     do m = 0,l
                        som = som + Qs_lmij(i,j,l+1,m+1,k)*plm(l,m,mu(n))
                     enddo
                     Qs_nij(i,j,n1,k) = Qs_nij(i,j,n1,k)   + 0.5D0*DBLE(2*l+1)*som     
               enddo
               n1=n1+1
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
                  Qs(i,j,k) = Qs(i,j,k) + 2.0D0*w(n)*Qs_nij(i,j,n1,k) !> end
                  n1=n1+1
               enddo
             enddo
           enddo
        enddo
        enddo 
        Qs_ij = Qs
        !print*, Qs_ij 
        !>---------------
        !do k = 1,ng
        !do j = 1,TNFM_y
        !   do i = 1,TNFM_x
        !      Qs_ij(i,j,k) = (0.5D0)*sum(SigS(fmmid2D(i,j),1,:,k)*phi_ij(i,j,:))
        !     
        !   enddo
        !enddo
        !enddo
    end subroutine Scattering_Source

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
