!BL1
subroutine gauleg(x1,x2,x,w,ngauss)
       implicit none
       integer(kind=4), intent(in) :: ngauss
       real(kind=8), intent(in) :: x1,x2
       real(kind=8), dimension(ngauss), intent(out) :: x,w
       real(kind=8), parameter :: EPSS = 3.0E-14
       !EPS is the relative precision.
       !Given the lower and upper limits of integration x1 and x2 ,
       !and given n ,   1rthis routine returns
       !arrays x(1:n) and w(1:n) of length n ,
       ! containing the abscissas and weights of the Gauss-
       !Legendre n-point quadrature formula.
       integer(kind=4) :: i,j,m
       real(kind=8) :: p1,p2,p3,pp,xl,xm,z,z1
       m=(ngauss+1)/2
       xm=0.50*(x2+x1)
       xl=0.50*(x2-x1)
       do  i= 1,m         !Loop over the desired roots.
       z=cos(3.1415926540*(i-.250)/(ngauss+.50))
       10 continue
       p1=1.0
       p2=0.0
          do  j=1,ngauss
          p3=p2
          p2=p1
          p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
          enddo 
       pp=ngauss*(z*p1-p2)/(z*z-1.0)
       z1=z
       z=z1-p1/pp
       if (abs(z-z1).gt.EPSS) goto 10
       x(i)=xm-xl*z
       x(ngauss+1-i)=xm+xl*z
       w(i)=2.0*xl/((1.0-z*z)*pp*pp)
       w(ngauss+1-i)=w(i)
       enddo 
       return     
end subroutine gauleg
!BL2
subroutine leg_poly(ngauss,order,mu,p) 
       implicit none
       integer(kind=4), intent(in) :: ngauss,order
       real(kind=8), dimension(ngauss), intent(in) :: mu
       real(kind=8), dimension(order,ngauss), intent(out) :: p
       real(kind=8) :: a1,a2
       integer(kind=4) :: l,k
       k = order
           if ( k == 1 ) then
              p(1,:) = 1.0
           elseif (k == 2 ) then
              p(1,:) = 1.0
              p(2,:) = mu(:)
           else
              p(1,:) = 1.0
              p(2,:) = mu(:)
              do l=2,order-1
                 a1 = (2*(float(l)-1.0)+1.0)/float(l)
                 a2 = (float(l)-1.0)/float(l)
                 p(l+1,:) = a1*mu(:)*p(l,:)-a2*p(l-1,:)
              enddo
           endif
end subroutine leg_poly
!BL3
subroutine fmm_id(totNFM,nregion,nfmesh,RegMat,fmmid)
       implicit none
       integer(kind=4), intent(in) :: totNFM,nregion
       integer(kind=4), dimension(nregion), intent(in) :: RegMat,nfmesh
       integer(kind=4), dimension(totNFM), intent(out) :: fmmid 
       integer(kind=4) :: m,i,j
       !-- fine mesh material id       
       m = 1
       do i = 1,size(nfmesh)
          do j=1,nfmesh(i)
             fmmid(m) = RegMat(i)
             m= m+1
          enddo
       enddo
end subroutine fmm_id
!BL4
subroutine Delta_f(nregion,totNFM,nfmesh,dcell,Delta)
       implicit none
       integer(kind=4), intent(in) :: nregion,totNFM
       integer(kind=4), dimension(nregion), intent(in) :: nfmesh
       real(kind=8), dimension(nregion), intent(in) :: dcell
       real(kind=8), dimension(totNFM), intent(out) :: Delta
       integer(kind=4) :: i,j,m
       m=1
       do i = 1,size(nfmesh)
          do j=1,nfmesh(i)
             Delta(m) = dcell(i)/nfmesh(i)
             m= m+1
          enddo
       enddo
end subroutine Delta_f
!BL5
subroutine Matrix_D(SigS,D,fmmid,ng,Nmat,order,totNFM,dim)
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,totNFM,dim,order
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
       real(kind=8), dimension(dim,dim,order), intent(out) :: D
       integer(kind=4) :: k0,k,k1,i,j 
       D = 0.0
       do k1 = 1,order
            k0 = 1
            do   k   = 1, ng
                 j   = 1
                 do  i  = k0,totNFM*k
                     D(i,i,k1) = SigS(fmmid(j),k1,k,k)
                     j = j + 1
                 end do
            k0 = totNFM + k0
            enddo
       enddo
end subroutine Matrix_D
!BL6
subroutine Matrix_L(SigS,L,fmmid,ng,Nmat,order,totNFM,dim)
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,totNFM,dim,order
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
       real(kind=8), dimension(dim ,dim,order), intent(out) :: L
       integer(kind=4) :: i,k0,k1,k2,k3,k4
                L(:,:,:) = 0.0
            do  k4  = 1,order
                k2  = 0
            do  k1  = 1,ng
                k0  = 1
                k3  = 0
                    do while (k0<k1) 
                             do  i  = 1,totNFM
                             L(i+k3,i+(k1-1)*totNFM,k4) = SigS(fmmid(i),k4,k0,k1) ! 1 --> 2 SigS(1,1,2)
                             enddo
                             k0 = k0 + 1
                             k3 = k3 + totNFM
                    enddo
            enddo
            enddo
end subroutine Matrix_L
!BL7
subroutine Matrix_U(SigS,U,fmmid,ng,Nmat,order,totNFM,dim)
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,totNFM,dim,order
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
       real(kind=8), dimension(dim,dim,order), intent(out) :: U
       integer(kind=4) :: i,k0,k1,k2,k3,k4
                U(:,:,:) = 0.0
            do  k4  = 1,order
                k2  = 0
            do  k1  = 1,ng
                k0  = 1
                k3  = 0
                    do while (k0<k1) 
                             do  i  = 1,totNFM
                             U(i+(k1-1)*totNFM,i+k3,k4) = SigS(fmmid(i),k4,k1,k0) ! 1 --> 2 SigS(1,1,2)
                             enddo
                             k0 = k0 + 1
                             k3 = k3 + totNFM
                    enddo
            enddo
            enddo
end subroutine Matrix_U
!BL8
subroutine Matrix_F(NusigF,Chi,F,fmmid,ng,Nmat,totNFM,dim)
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,totNFM,dim
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=8), dimension(Nmat,ng), intent(in) :: NusigF,Chi
       real(kind=8), dimension(dim ,dim), intent(out) :: F
       integer(kind=4) :: i,j,k0,k1,k2,k3,k
            F(:,:) = 0.0
            k0 = 1
            do   k   = 1, ng
                 j   = 1
                 do  i  = k0,totNFM*k
                     F(i,i) = Chi(fmmid(j),k)*NusigF(fmmid(j),k)
                     j = j + 1
                 end do
            k0 = totNFM + k0
            enddo
              
            k2  = 0
            do  k1  = 1,ng
                k0  = 1
                k3  = 0
                    do while (k0<k1) 
                             do  i  = 1,totNFM
                             F(i+k3,i+(k1-1)*totNFM) =  Chi(fmmid(i),k1)*NusigF(fmmid(i),k0)
                             F(i+(k1-1)*totNFM,i+k3) =  Chi(fmmid(i),k0)*NusigF(fmmid(i),k1)
                             enddo
                             k0 = k0 + 1
                             k3 = k3 + totNFM
                    enddo
            enddo
end subroutine Matrix_F
!BL9
subroutine flux_guess(dim,ng,Nmat,nregion,ngauss,order,totNFM,NusigF,dcell,p,flux_ni,flux_li)
       implicit none
       integer(kind=4), intent(in) :: dim, ng,Nmat,nregion,ngauss,order,totNFM
       real(kind=8), dimension(Nmat,ng), intent(in) :: NusigF
       real(kind=8), dimension(order,ngauss), intent(in) :: p
       real(kind=8), dimension(totNFM), intent(in) :: dcell
       real(kind=8), dimension(dim,ngauss*ng), intent(out) :: flux_ni
       real(kind=8), dimension(dim,order), intent(out) :: flux_li
       integer(kind=4) :: i,k1,n,ll,m,k0,k2
       real(kind=8), dimension(ng*nregion) :: a10,a11
       i = 1
       do k1=1,ng
          do k2=1,nregion
             a10(i) = NusigF(k2,k1)
             a11(i) = dcell(k2)
             i = i + 1
          enddo
       enddo

       flux_li(:,:)  = 1.0/dot_product(a10,a11)
       k2=1
       k0=1
       do k1 = 1,ng
          do i =k0,totNFM*k1
              do n = 1,ngauss
                ll = 0
                do m = 1,order
                   flux_ni(i,n+k2-1) = flux_ni(i,n+k2-1) + 0.5*(2.*float(ll) +1.0)*flux_li(i,m)*p(m,n) 
                   ll = ll+1
                enddo
              enddo   
          enddo
             k2 = k2 + ngauss
             k0 = k0 + totNFM
       enddo
       flux_ni = 1
end subroutine flux_guess
!BL10
subroutine Fission_Source(ng,dim,ngauss,order,totNFM,F,flux_li,p,k_eff,FQ_ni)
       implicit none
       integer(kind=4), intent(in) :: dim,ngauss,ng,totNFM,order
       real(kind=8), intent(in) :: k_eff 
       real(kind=8), dimension(dim ,dim), intent(in) :: F
       real(kind=8), dimension(order,ngauss), intent(in) :: p
       real(kind=8), dimension(dim,order), intent(in) :: flux_li
       real(kind=8), dimension(dim,ngauss*ng), intent(out) :: FQ_ni
       real(kind=8), dimension(dim,order) :: Q_li,phi_li
       integer(kind=4) :: i,k0,k1,k2,ll,n,m
       k0=1;k2=1
       phi_li = flux_li
       FQ_ni(:,:) = 0.0
       if (order >= 2) then  
          do i = 2, order
             phi_li(:,i) = 0.0
          enddo
       endif

       Q_li(:,:) = matmul(F(:,:),phi_li(:,:))

       do k1 = 1,ng
          do i =k0,totNFM*k1
              do n = 1,ngauss
                ll = 0
                do m = 1,order
                    FQ_ni(i,n+k2-1) = FQ_ni(i,n+k2-1) + 0.5*(2.*float(ll)+1.)*Q_li(i,m)*p(m,n)/k_eff
                   ll = ll+1
                enddo
              enddo   
          enddo
             k2 = k2 + ngauss
             k0 = k0 + totNFM
       enddo
end subroutine Fission_Source
!BL11
subroutine Scattering_Source(ng,dim,ngauss,order,totNFM,D,U,L,flux_li,p,SQ_ni)
       implicit none
       integer(kind=4), intent(in) :: dim,ngauss,ng,totNFM,order
       real(kind=8), dimension(dim,dim,order), intent(in) :: L,U,D
       real(kind=8), dimension(order,ngauss), intent(in) :: p
       real(kind=8), dimension(dim,order), intent(in) :: flux_li
       real(kind=8), dimension(dim,ngauss*ng), intent(out) :: SQ_ni
       real(kind=8), dimension(dim,order) :: Q_li
       integer(kind=4) :: i,k0,k1,k2,ii,ll,n,m
       SQ_ni(:,:) = 0.0

       do ii = 1,order
          Q_li(:,ii) = matmul(D(:,:,ii),flux_li(:,ii)) + &
                       matmul(L(:,:,ii),flux_li(:,ii)) + &
                       matmul(U(:,:,ii),flux_li(:,ii))
       enddo

       k2=1
       k0=1
       do k1 = 1,ng
          do i =k0,totNFM*k1
              do n = 1,ngauss
                ll = 0
                do m = 1,order
                   SQ_ni(i,n+k2-1) = SQ_ni(i,n+k2-1) + 0.5*(2.*float(ll) +1.)*Q_li(i,m)*p(m,n) 
                   ll = ll+1
                enddo
              enddo   
          enddo
             k2 = k2 + ngauss
             k0 = k0 + totNFM
       enddo
end subroutine Scattering_Source
!BL12
subroutine Total_Source(ng,dim,ngauss,FQ_ni,SQ_ni,Q_ni)
       implicit none
       integer(kind=4), intent(in) :: dim,ngauss,ng
       real(kind=8), dimension(dim,ngauss*ng), intent(in) :: FQ_ni,SQ_ni
       real(kind=8), dimension(dim,ngauss*ng), intent(out) :: Q_ni
       
       Q_ni =  SQ_ni + FQ_ni
end subroutine Total_Source
!BL13
subroutine Matrix_AB(ng,Nmat,dim,totNFM,ngauss,mu,fmmid,SigT,Delta,Q_ni,A,B)
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,dim,totNFM,ngauss
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=8), dimension(totNFM), intent(in) :: Delta
       real(kind=8), dimension(ngauss), intent(in) :: mu
       real(kind=8), dimension(Nmat,ng), intent(in) :: SigT
       real(kind=8), dimension(dim,ngauss*ng), intent(in) :: Q_ni
       real(kind=8), dimension(dim,ngauss*ng), intent(out) :: A,B
       integer(kind=4) :: i,n,k0,k1,k2,j
       k0=1;k2=1
     
       do k1=1,ng
                 j = totNFM
             do  i = totNFM*k1,k0,-1 ! right-to-left 
                 do  n = 1,(ngauss/2)
                     A(i,n + k2 -1) = (SigT(fmmid(j),k1)*Delta(j))/abs(mu(n))
                     B(i,n + k2 -1) = Q_ni(i,n + k2 -1)/SigT(fmmid(j),k1)
                 enddo
                 j = j-1
             enddo

                     
                 j = 1
             do  i = k0,totNFM*k1    ! left-to-rights
                 do  n = ngauss,(ngauss/2)+1,-1
                     A(i,n + k2 -1) = (SigT(fmmid(j),k1)*Delta(j))/abs(mu(n))
                     B(i,n + k2 -1) = Q_ni(i,n + k2 -1)/SigT(fmmid(j),k1)
                 enddo
                     j = j + 1
             enddo

             k0 = k0 + totNFM
             k2 = k2 + ngauss
       enddo
end subroutine Matrix_AB
!BL14
subroutine current_f(BC,scheme,A,B,ng,ngauss,totNFM,dim,curr)
       implicit none
       CHARACTER(50), intent(in) :: BC, scheme
       integer(kind=4), intent(in) :: ng,ngauss,totNFM,dim
       real(kind=8), dimension(dim,ngauss*ng), intent(in) :: A,B
       real(kind=8), dimension((totNFM+1)*ng,ngauss*ng), intent(out) :: curr
       real(kind=8), dimension((totNFM+1)*ng,ngauss*ng) ::curra,currb
       real(kind=8) :: al
       integer(kind=4) :: i,j,k1,k0,k2,m,n
       curr = 0.0
       k0=1;k2=1;n = 1
       !========================\\vacuum//==============================
       if  ( BC == 'Vacuum' ) then 
           do k1 = 1,ng
                        m = k1*totNFM
                do  i = k1*totNFM+k1-1,k0,-1 ! right-to-left
                    do  j=1,(ngauss/2)
                        if (scheme == 'DD0') then
                        al = (2.-A(m,j+k2-1))/(2.+A(m,j+k2-1))
                        elseif (scheme == 'DD1') then
                        al = (12.-6*A(m,j+k2-1)+A(m,j+k2-1)**2)/(12.+6*A(m,j+k2-1)+A(m,j+k2-1)**2)
                        else
                        al = exp(-A(m,j+k2-1))
                        endif
                        curr(i,j+k2-1) = curr(i+1,j+k2-1)*al + B(m,j+k2-1)*(1-al)
                    enddo
                        m = m - 1
                enddo    

                        n = m + 1 
                do  i = k0,k1*totNFM+k1-1   ! left-to-right
                    do  j=ngauss,(ngauss/2)+1,-1
                        if (scheme == 'DD0') then
                        al = (2.-A(n,j+k2-1))/(2.+A(n,j+k2-1))
                        elseif (scheme == 'DD1') then
                        al = (12.-6*A(n,j+k2-1)+A(n,j+k2-1)**2)/(12.+6*A(n,j+k2-1)+A(n,j+k2-1)**2)
                        else
                        al = exp(-A(n,j+k2-1))
                        endif
                        curr(i+1,j+k2-1) = curr(i,j+k2-1)*al + B(n,j+k2-1)*(1-al)
                    enddo
                        n = n + 1
                enddo
                        k0 = k0 + (totNFM + 1)
                        k2 = k2 + ngauss
           enddo
       !========================\\reflective//==============================

        elseif  ( BC == 'Reflective' ) then 
            do k1 = 1,ng
                        m = k1*totNFM
                do  i = k1*totNFM+k1-1,k0,-1 ! right-to-left 
                    do  j=1,(ngauss/2)
                        if (scheme == 'DD0') then
                        al = (2.-A(m,j+k2-1))/(2.+A(m,j+k2-1))
                        elseif (scheme == 'DD1') then
                        al = (12.-6*A(m,j+k2-1)+A(m,j+k2-1)**2)/(12.+6*A(m,j+k2-1)+A(m,j+k2-1)**2)
                        else
                        al = exp(-A(m,j+k2-1))
                        endif
                        curra(k1*totNFM+k1,j + k2 -1) = 0.0
                        currb(k1*totNFM+k1,j + k2 -1) = 1.0
                        curra(i,j+k2-1) = curra(i+1,j+k2-1)*al + B(m,j+k2-1)*(1-al)
                        currb(i,j+k2-1) = currb(i+1,j+k2-1)*al + B(m,j+k2-1)*(1-al)
                        !   On checrche le courant dans la cellule 
                        curr(i+1,j+k2-1)= curra(i,j+k2-1)/ (1.0 + (curra(i,j+k2-1)-currb(i,j+k2-1)))
                    enddo
                        m = m - 1
                enddo  

                        n = m + 1
                do  i = k0,k1*totNFM+k1-1    ! left-to-right
                    do  j=ngauss,(ngauss/2)+1,-1
                        curra(k0,j+k2-1) = 0.0
                        currb(k0,j+k2-1) = 1.0
                        if (scheme == 'DD0') then
                        al = (2.-A(n,j+k2-1))/(2.+A(n,j+k2-1))
                        elseif (scheme == 'DD1') then
                        al = (12.-6*A(n,j+k2-1)+A(n,j+k2-1)**2)/(12.+6*A(n,j+k2-1)+A(n,j+k2-1)**2)
                        else
                        al = exp(-A(n,j+k2-1))
                        endif
                        curra(i+1,j+k2-1) = curra(i,j+k2-1)*al + B(n,j+k2-1)*(1-al)
                        currb(i+1,j+k2-1) = currb(i,j+k2-1)*al + B(n,j+k2-1)*(1-al)
                        !   Conditions au limite
                        curr(i,j+k2-1)=curra(i+1,j+k2-1)/(1.0 + (curra(i+1,j+k2-1)-currb(i+1,j+k2-1)))  
                    enddo
                        n = n + 1
                enddo
                        k0 = k0 + (totNFM + 1)
                        k2 = k2 + ngauss
           enddo

           k0=1;k2=1

           do k1 = 1,ng
                        m = k1*totNFM
                do  i = k1*totNFM+k1-1,k0,-1 ! right-to-left 
                    do  j=1,(ngauss/2)
                        if (scheme == 'DD0') then
                        al = (2.-A(m,j+k2-1))/(2.+A(m,j+k2-1))
                        elseif (scheme == 'DD1') then
                        al = (12.-6*A(m,j+k2-1)+A(m,j+k2-1)**2)/(12.+6*A(m,j+k2-1)+A(m,j+k2-1)**2)
                        else
                        al = exp(-A(m,j+k2-1))
                        endif
                        curr(i,j+k2-1) = curr(i+1,j+k2-1)*al + B(m,j+k2-1)*(1-al)
                    enddo
                        m = m - 1
                enddo  
                        n = m + 1
                do  i = k0,k1*totNFM+k1-1    ! left-to-right
                    do  j=ngauss,(ngauss/2)+1,-1
                        if (scheme == 'DD0') then
                        al = (2.-A(n,j+k2-1))/(2.+A(n,j+k2-1))
                        elseif (scheme == 'DD1') then
                        al = (12.-6*A(n,j+k2-1)+A(n,j+k2-1)**2)/(12.+6*A(n,j+k2-1)+A(n,j+k2-1)**2)
                        else
                        al = exp(-A(n,j+k2-1))
                        endif
                        curr(i+1,j+k2-1) = curr(i,j+k2-1)*al + B(n,j+k2-1)*(1-al)
                    enddo
                        n = n + 1
                enddo
                        k0 = k0 + (totNFM + 1)
                        k2 = k2 + ngauss
           enddo
           !========================\\reflective_vacuum//==============================
       elseif  ( BC == 'Vacuum Reflective' ) then 
             do k1 = 1,ng
                    n = (k1-1)*totNFM + 1
                do  i = k0,k1*totNFM+k1-1    ! left-to-right
                    do  j=(ngauss/2+1),ngauss
                        curr(k0,j + k2 -1) = 0.0
                        if (scheme == 'DD0') then
                        al = (2.-A(n,j+k2-1))/(2.+A(n,j+k2-1))
                        elseif (scheme == 'DD1') then
                        al = (12.-6.*A(n,j+k2-1)+A(n,j+k2-1)**2)/(12.+6.*A(n,j+k2-1)+A(n,j+k2-1)**2)
                        else
                        al = exp(-A(n,j+k2-1))
                        endif
                        curr(i+1,j+k2-1) =  curr(i,j+k2-1)*al + B(n,j+k2-1)*(1-al)             
                    enddo
                    n = n + 1
                enddo
   
                    m = n - 1
                do  i = k1*totNFM+k1-1,k0,-1 ! right-to-left
                    do  j=1,(ngauss/2)
                        if (scheme == 'DD0') then
                        al = (2.-A(m,j+k2-1))/(2.+A(m,j+k2-1))
                        elseif (scheme == 'DD1') then
                        al = (12.-6*A(m,j+k2-1)+A(m,j+k2-1)**2)/(12.+6*A(m,j+k2-1)+A(m,j+k2-1)**2)
                        else
                        al = exp(-A(m,j+k2-1))
                        endif
                        curr(k1*totNFM+k1,j+k2 -1) = curr(k1*totNFM+k1,ngauss-j+k2)
                        curr(i,j+k2-1) = curr(i+1,j+k2-1)*al + B(m,j+k2-1)*(1-al)
                    enddo
                    m = m - 1
                enddo 
                k0 = k0 + (totNFM + 1)
                k2 = k2 + ngauss
             enddo
       elseif  ( BC == 'Reflective Vacuum' ) then
             do k1 = 1,ng
                    m = k1*totNFM
                do  i = k1*totNFM+k1-1,k0,-1 ! right-to-left 
                    do  j=(ngauss/2),1,-1
                        if (scheme == 'DD0') then
                        al = (2.-A(m,j+k2-1))/(2.+A(m,j+k2-1))
                        elseif (scheme == 'DD1') then
                        al = (12.-6.*A(m,j+k2-1)+A(m,j+k2-1)*A(m,j+k2-1))/(12.&
                             +6.*A(m,j+k2-1)+A(m,j+k2-1)*A(m,j+k2-1))
                        else
                        al = exp(-A(m,j+k2-1))
                        endif
                        curr(k1*totNFM+k1,j + k2 -1) = 0.0
                        curr(i,j+k2-1) = curr(i+1,j+k2-1)*al + B(m,j+k2-1)*(1-al)
                    enddo
                    m = m -1
                enddo
                    n = m + 1
                do  i = k0,k1*totNFM+k1-1    ! left-to-right
                    do  j=(ngauss/2+1),ngauss
                        if (scheme == 'DD0') then
                        al = (2.-A(n,j+k2-1))/(2.+A(n,j+k2-1))
                        elseif (scheme == 'DD1') then
                        al = (12.-6.*A(n,j+k2-1)+A(n,j+k2-1)*A(n,j+k2-1))/(12.&
                             +6.*A(n,j+k2-1)+A(n,j+k2-1)*A(n,j+k2-1))
                        else
                        al = exp(-A(n,j+k2-1))
                        endif
                        curr(k0,j + k2 -1) = curr(k0,ngauss-j+k2) 
                        curr(i+1,j+k2-1) = curr(i,j+k2-1)*al + B(n,j+k2-1)*(1-al)
                    enddo
                    n = n + 1
                enddo
                k0 = k0 + (totNFM + 1)
                k2 = k2 + ngauss
             enddo
       endif
end subroutine current_f
!BL15
subroutine flux_moy(A,B,ng,ngauss,totNFM,dim,curr,flux_ni)
       implicit none
       integer(kind=4), intent(in) :: ng,ngauss,totNFM,dim
       real(kind=8), dimension(dim,ngauss*ng), intent(in) :: A,B
       real(kind=8), dimension((totNFM+1)*ng,ngauss*ng), intent(in) :: curr
       real(kind=8), dimension(dim,ngauss*ng), intent(out) :: flux_ni
       integer(kind=4) :: i,j,k1,k0,k2,m,n,l       
                        flux_ni = 0.0
                        k0 = 1
                        k2 = 1
                        l = 0
            do k1 = 1,ng
                        m = k1*totNFM
                do  i = totNFM*k1,k0,-1 ! right-to-left 
                    do  j = 1,(ngauss/2)
                        flux_ni(i,j+k2-1) = (curr(i+l+1,j+k2-1)-curr(i+l,j+k2-1))/A(m,j+k2-1) + B(m,j+k2-1)  
                    enddo
                        m = m - 1
                enddo  
                        n = m + 1 
                do  i = k0,k1*totNFM    ! left-to-right
                    do  j= ngauss,(ngauss/2)+1,-1
                        flux_ni(i,j+k2-1) =  (curr(i+l,j+k2-1)-curr(i+l+1,j+k2-1))/A(n,j+k2-1) + B(n,j+k2-1) 
                    enddo
                        n = n + 1
                enddo
                        k0 = k0 + totNFM 
                        k2 = k2 + ngauss
                        l = l + 1
            enddo
end subroutine flux_moy
!BL16
subroutine Outer_Iteration(ng,dim,Max_it,totNFM,ngauss,order,Nmat,it,inter,eps,wt,mu,D,F,U,L,p,BC,scheme,fmmid,&
                        SigT,flux_ni,flux_li,Delta,k_eff,phi)
       implicit none
       integer(kind=4), intent(in) :: ng,dim,totNFM,ngauss,order,Nmat,Max_it
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid 
       real(kind=8), dimension(ngauss), intent(in) :: wt,mu
       real(kind=8), dimension(totNFM), intent(in) :: Delta
       real(kind=8), dimension(dim ,dim), intent(in) :: F
       real(kind=8), dimension(order,ngauss), intent(in) :: p
       real(kind=8), dimension(Nmat,ng), intent(in) :: SigT
       real(kind=8), dimension(dim,dim,order), intent(in) :: D,U,L
       real(kind=8), intent(in) :: eps
       real(kind=8), dimension(dim,ngauss*ng), intent(inout) :: flux_ni
       real(kind=8), dimension(dim,order), intent(inout) :: flux_li
       real(kind=8), dimension(dim), intent(out) :: phi
       real(kind=8), intent(out) :: k_eff
       integer(kind=4), intent(out) :: it,inter
!      variables locale
       real(kind=8), dimension(dim,ngauss*ng) :: A,B
       real(kind=8), dimension(dim,ngauss*ng) :: Q_ni,SQ_ni,flux0,flux00,flux,FQ_ni
       real(kind=8), dimension((totNFM+1)*ng,ngauss*ng) :: curr
       real(kind=8), dimension(dim,order) :: flux_li0
       real(kind=8), dimension(ng) :: moy 
       real(kind=8) :: err_k_eff, err_phi, k_eff0, eval1, eval2,Del,Sig,muu,err_flux
       integer(kind=4) :: i,n,k0,k2,k1
       CHARACTER(50), intent(in) :: BC, scheme
!      convergence parameters
       inter=0
       err_k_eff = 1.0
       err_phi = 1.0
       err_flux = 1.5
       it = 0 
       k_eff = 1.0
!      the positive flux condition
       Del = minval(Delta)
       Sig = minval(SigT)
       muu = minval(mu) 
       
       if (Del*Sig > 2*abs(muu)) then
          print*,'Failed the positive flux condition.'
          stop
          endif 
       
       flux_li0 = flux_li
       flux0    = abs(flux_ni)  
       flux00   = abs(flux_ni)
            
       do while ( err_k_eff >= eps .and. err_phi >= eps )
                 flux0 = abs(flux_ni)
                 call Fission_Source(ng,dim,ngauss,order,totNFM,F,flux_li,p,k_eff,FQ_ni)
                 if (it >= Max_it) then
                 print*,'Failed to converge.'
                 stop
                 endif 
                 flux = 0.0
                 k_eff0 = k_eff
                 err_flux = 1.5
! Staring intern Iteration
! ==============================================================================
             do  while (  err_flux >= 0.0001 )
                 
                 call Scattering_Source(ng,dim,ngauss,order,totNFM,D,U,L,flux_li,p,SQ_ni)
                 call Total_Source(ng,dim,ngauss,FQ_ni,SQ_ni,Q_ni)
                 call Matrix_AB(ng,Nmat,dim,totNFM,ngauss,mu,fmmid,SigT,Delta,Q_ni,A,B)
                 call current_f(BC,scheme,A,B,ng,ngauss,totNFM,dim,curr)
                 call flux_moy(A,B,ng,ngauss,totNFM,dim,curr,flux_ni)

                 flux_li = 0.0
                 k0 = 1
                 k2 = 1
                 do k1 = 1,ng
                    do i  = k1*totNFM,k0,-1
                       do n = 1,(ngauss/2)
                          flux_li(i,:) = flux_li(i,:) + wt(n)*abs(flux_ni(i,n+k2-1))*p(:,n)
                       enddo
                    enddo 

                    do i = k0,totNFM*k1
                       do n = ngauss,(ngauss/2)+1,-1
                           flux_li(i,:) = flux_li(i,:) + wt(n)*abs(flux_ni(i,n+k2-1))*p(:,n)
                       enddo
                    enddo
                       k0 = k0 + totNFM
                       k2 = k2 + ngauss
                 enddo

                 flux =  flux + abs(flux_ni)  
                 err_flux = maxval(abs(flux-flux0)/flux0)
                 flux0 = flux 
                 inter = inter + 1 
             enddo
! ending intern Iteration 
!==============================================================================
             phi = 0.0
                 k0 = 1
                 k2 = 1
                 do k1 = 1,ng
                    do i  = k1*totNFM,k0,-1
                       do n = 1,(ngauss/2)
                          phi(i) = phi(i) + wt(n)*abs(flux_ni(i,n+k2-1))
                       enddo
                    enddo 

                    do i = k0,totNFM*k1
                       do n = ngauss,(ngauss/2)+1,-1
                           phi(i) =  phi(i) + wt(n)*abs(flux_ni(i,n+k2-1))
                       enddo
                    enddo
                       k0 = k0 + totNFM
                       k2 = k2 + ngauss
                 enddo

             err_phi =  maxval(abs(abs(flux_ni)-flux00)/flux00)
             it = it + 1     
             eval1=sum(matmul(F, flux_li(:,1)))
             eval2=sum(matmul(F, flux_li0(:,1)))
             k_eff = k_eff*(eval1/eval2)
             err_k_eff =  abs(k_eff-k_eff0)/k_eff
             write(*,2000)it,k_eff,err_k_eff  
             flux_li0 = flux_li
             flux00   = abs(flux_ni)
       end do  
       moy = 0.0
       k1 = 1
       do k0 = 1,ng
          do i=k1,totNFM*k0
             moy(k0) = moy(k0) + phi(i)
          end do
          do i=k1,totNFM*k0
             phi(i) = phi(i)/(moy(k0)/totNFM)
          end do
          k1 = k1 + totNFM 
       enddo
       2000 format("Iteration",i4,":",1x,"===>",1x,"keff =",F9.6,1x,"===>",1x,"res =",e10.3)
end subroutine Outer_Iteration
!BL17
subroutine Output(start,BC,tm,k_eff,SigT,NusigF,SigS,Chi,mu,wt,dcell,phi,eps,totNFM,dim,&
                  ng,Nmat,order,nregion,ngauss,it1,it2)
        implicit none
        integer(kind=4), intent(in) :: ng,dim,totNFM,Nmat,order,nregion,ngauss,it1,it2
        real(kind=8), dimension(Nmat,ng), intent(in) :: SigT,NusigF,Chi
        real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
        real(kind=8), dimension(dim), intent(in) :: phi
        real(kind=8), dimension(nregion), intent(in) :: dcell
        real(kind=8), dimension(ngauss), intent(in) :: mu,wt
        CHARACTER(50), intent(in) :: start,BC,tm
        real(kind=8), intent(in) :: eps,k_eff
        integer(kind=4) :: i,j
        open (100,file='app/Output/OUTPUT_MOC.TXT')
        write (100, FMT=* ) '********************************************************************************'
        write (100, FMT=* ) 'ERSN, UNIVERSITY ABDELMALEK ESSAADI FACULTY OF SCIENCES - TETOUAN, MOROCCO'
        write (100, FMT=* ) 'CODE  DEVELOPED  BY  MOHAMED  LAHDOUR,  PHD  STUDENT'
        write (100, FMT=* ) 'NTP-ERSN:        MOC  MOETHOD OF CHARACTERISTICS'
        write (100, FMT=* ) 'VERSION NUMBER:  1.2'
        write (100, FMT=* ) 'VERSION DATE:    8  OTOBER  2018'
        write (100,3010) 'RAN ON:          ', start,'(H:M:S)'
        write (100, FMT=* ) '********************************************************************************'
        write (100, FMT=* ) '           ----------------------------------------------------------' 
        write (100, FMT=* ) '                     INPUT  PARAMETER - VALUES  FROM  INPUT'              
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'ENERGY GROUP NUMBER:                   ',ng
        write (100, FMT=* ) 'REGIONS NUMBER:                        ',nregion
        write (100, FMT=* ) 'MATERIALS NUMBER:                      ',Nmat
        write (100,3040)    'SIZE FOR EACH MATERIAL PER [CM]:       ',dcell       
        write (100, FMT=* ) 'DISCRETIZATIONS ANGULAR:               ',ngauss
        write (100, FMT=* ) 'L-ORDER LEGENDRE POLYNOMIAL:           ',order-1
        write (100, FMT=* ) 'TOTAL NUMBER OF FINE MESHES:           ',totNFM
        write (100,3050)    'CONVERGENCE CRITERION of KEFF AND FLUX:',eps
        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '                      CALCULATION  RUN-TIME  PARAMETERS  MOC' 
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'GAUSS  LEGENDRE  QUADRATURE  POINTS  AND  WEIGHTS: '
        write (100, FMT=* ) ''
        write (100, FMT=* ) '      N. GAUSS ','         POINTS    ','     WEIGHTS '
        write (100, FMT=* ) ''
        do i=1,ngauss
        write(100,3060) i,mu(i),wt(i)
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
        write (100, FMT=* ) '                             SCALAR  FLUX  SOLUTION' 
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'FLUXES  PER  MESH  PER  ENERGY  GROUP:'  
        write (100, FMT=* ) '' 
        write (100,3000)'       M E S H ', ('     G R O U P',i,i=1,ng)
        write (100, FMT=* ) ''
        do i=1,totNFM
        write(100,2000) i,(phi(i+j), j=0,dim-1,totNFM)
        enddo
        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '             OUTPUT  PARAMETER - SOLUTION  TO  TRANSPORT  EQUATION' 
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        if  ( BC == 'reflective' ) then
        write (100,3090)    'K-INF                    =',k_eff
        else
        write (100,3090)    'K-EFF                    =',k_eff
        endif
        write (100,3020)    'N. OUTER ITERATIONS      =',it1
        write (100,3020)    'TOTAL INNER ITERATIONS   =',it2
        write (100,4000)    'TOTAL EXECUTION TIME     =',tm,'(H:M:S)'
        write (100, FMT=* ) ''
        write (100, FMT=* ) '********************************************************************************'
        2000 format(1x,1p,i11,5x,200e16.5) 
        3000 format(1x,A14,2x,300(A14,i2))  
        3010 format(1x,A17,A22,A10)
        3020 format(1x,A26,4x,i10)
        3040 format(1x,A33,2x,200F10.5)
        3050 format(1x,1p,A41,4x,e8.1)
        3060 format(1x,1p,i11,5x,e16.5,e16.5)
        3070 format(1x,A18,i4)
        3080 format(1x,1p,i11,5x,e16.5,e16.5,e16.5,e16.5,e16.5)
        3090 format(1x,A26,6x,f8.6)
        4000 format(1x,A26,4x,A10,A10)
        close(100)
end subroutine Output
!BL18
subroutine plot_flux(Delta,flux,totNFM,dim)
        implicit none
        integer(kind=4), intent(in) :: dim,totNFM
        real(kind=8), dimension(dim), intent(in) :: flux
        real(kind=8), dimension(totNFM), intent(in) :: Delta
        real(kind=8) :: som
        integer(kind=4) :: i,j
        open (10,file='app/Output/flux_moc.h')
        som = Delta(1)
        do i=1,totNFM
        write(10,*) som,(flux(i+j),j=0,dim-1,totNFM) !/flux(1+j)
        som = som + Delta(i)
        enddo
        close(10)
        !2000 format(1x,100p6e12.4)
end subroutine plot_flux
!BL19
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
       print*, 'Copyright| 2015-2019 Radiation and Nuclear system Laboratory'
       print*, '         | University Abdelmalk Essaadi, FS Tetouan, Morocco'
       print*, 'Source   | FORTRAN90 version                                ' 
       print*, 'GUI      | PyQt5                                            ' 
       print*, 'Method   | The Method of Characteristics (MOC)              '  
       print*, 'Dimension| One dimensions (1D)                              '
       print*, 'Geometry | Slab                                             '  
       print*, '____________________________________________________________'
       print*, ''
    end subroutine title1
!BL20
    subroutine title2()
       write(*,FMT=*)'**********************************************************'
       write(*,FMT=*)'                         Finished                         '                             
       write(*,FMT=*)'**********************************************************'  
    end subroutine title2
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

