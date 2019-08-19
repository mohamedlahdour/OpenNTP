!BL1
subroutine gauleg(x1,x2,mu,w,ngauss)
       implicit none
       integer(kind=4), intent(in) :: ngauss
       real(kind=8), intent(in) :: x1,x2
       real(kind=8), dimension(ngauss), intent(out) :: mu,w
       real(kind=8), dimension(ngauss) :: x
       real(kind=8), parameter :: EPSS = 3.0d-14
       !EPS is the relative precision.
       !Given the lower and upper limits of integration x1 and x2 ,
       !and given n , this routine returns
       !arrays x(1:n) and w(1:n) of length n ,
       ! containing the abscissas and weights of the Gauss-
       !Legendre n-point quadrature formula.
       integer(kind=4) :: i,j,m
       real(kind=8) :: p1,p2,p3,pp,xl,xm,z,z1
       m=(ngauss+1)/2
       xm=0.5d0*(x2+x1)
       xl=0.5d0*(x2-x1)
       do  i=1,m         !Loop over the desired roots.
       z=cos(3.141592654d0*(i-.25d0)/(ngauss+.5d0))
       10 continue
       p1=1.d0
       p2=0.d0
          do  j=1,ngauss
          p3=p2
          p2=p1
          p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
          enddo 
       pp=ngauss*(z*p1-p2)/(z*z-1.d0)
       z1=z
       z=z1-p1/pp
       if (abs(z-z1).gt.EPSS) goto 10
       x(i)=xm-xl*z
       x(ngauss+1-i)=xm+xl*z
       w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
       w(ngauss+1-i)=w(i)
       enddo 
       mu=-x
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
                 a1 = (2*(real(l)-1.0)+1.0)/real(l)
                 a2 = (real(l)-1.0)/real(l)
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
subroutine Matrix_AB(ng,Nmat,dim,totNFM,ngauss,mu,fmmid,SigT,Delta,A,B)
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,dim,totNFM,ngauss
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=8), dimension(totNFM), intent(in) :: Delta
       real(kind=8), dimension(ngauss), intent(in) :: mu
       real(kind=8), dimension(Nmat,ng), intent(in) :: SigT
       real(kind=8), dimension(dim,ngauss*ng), intent(out) :: A,B
       integer(kind=4) :: i,n,k0,k1,k2,j
       k0=1;k2=1 
       do k1=1,ng
                 j = 1
             do  i = k0,totNFM*k1    ! sweeping to the right
                 do  n = 1,ngauss/2
                     A(i,n + k2 -1) = 2.d0*mu(n)/Delta(j)+SigT(fmmid(j),k1)
                     B(i,n + k2 -1) = -2.0d0*mu(n)/Delta(j)
                     !print*,A(i,n + k2 -1),i,n + k2 -1
                     !print*,B(i,n + k2 -1)
                 enddo
                     j = j + 1
             enddo
                     
                 j = totNFM
             do  i = k0,totNFM*k1    ! sweeping to the left
                 do  n = (ngauss/2)+1,ngauss
                     !print*, mu(n)
                     A(i,n + k2 -1) = -2.d0*mu(n)/Delta(j)+SigT(fmmid(j),k1) 
                     B(i,n + k2 -1) =  2.d0*mu(n)/Delta(j)
                 enddo
                     j = j - 1
             enddo
             k0 = k0 + totNFM
             k2 = k2 + ngauss
       enddo
end subroutine Matrix_AB
!BL10
subroutine Matrix_MP(ng,Nmat,dim,totNFM,ngauss,A,B,MP,MM)
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,dim,totNFM,ngauss
       real(kind=8), dimension(dim,ngauss*ng), intent(in) :: A,B
       real(kind=8), dimension(dim,dim,ngauss/2), intent(out) :: MP,MM
       integer(kind=4) :: i,j,n,k0,k1,k2,k3,k
       k0=1;k2=1;k3=2
       do k1=1,ng
             do  n = 1,ngauss/2          ! sweeping to the right 
                 do  i = k0,totNFM*k1                   
                     MP(i,i,n) = A(i,n)           
                 enddo
                 do j=k0,totNFM*k1-1
                    k=1
                    do i = j+1,totNFM*k1
                       If (Mod(j,2).Eq.0) Then
                       MP(i,k,n) =  -2.0d0*B(i,n) 
                       k=k+1 
                       else
                       MP(i,k,n) =  2.0d0*B(i,n) 
                       k=k+1
                       endif
                    enddo
                 enddo
             enddo   
             do  n = (ngauss/2)+1,ngauss ! sweeping to the left
                 do  i = k0,totNFM*k1                   
                     MM(i,i,n-ngauss/2) = A(i,n)           
                 enddo

                 do j=k0,totNFM*k1-1
                     k=1
                    do i = j+1,totNFM*k1
                       If (Mod(j,2).Eq.0) Then
                       MM(i,k,n-ngauss/2) =  -2.0d0*B(i,n)
                       k=k+1 
                       else
                       MM(i,k,n-ngauss/2) =  2.0d0*B(i,n)
                       k=k+1 
                       endif
                    enddo
                 enddo
             enddo 
  
             k0 = k0 + totNFM
             k2 = k2 + ngauss
       enddo
       !write(*,'(5f8.4)') transpose(MP(:,:,1))
       !write(*,'(5f8.4)') transpose(MM(:,:,2))
end subroutine Matrix_MP
!BL11
subroutine Matrix_A(ng,Nmat,dim,totNFM,ngauss,MP,MM,A)
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,dim,totNFM,ngauss
       real(kind=8), dimension(dim,dim,ngauss/2), intent(in) :: MP,MM
       real(kind=8), dimension(dim*ngauss,dim*ngauss), intent(out) :: A
       real(kind=8), dimension(dim*ngauss,dim*ngauss) :: B
       integer(kind=4) :: i,j,n
       open (100,file='app/Output/matrixa_sn.h')
       n=1; B=0.0D0; A=0.0D0
       do i=1,dim*ngauss/2,totNFM
          B(i:i-1+totNFM,i:i-1+totNFM)=MP(1:totNFM,1:totNFM,n)
          n=n+1
       enddo
         n=1
       do i=dim*ngauss/2+1,dim*ngauss,totNFM
          B(i:i-1+totNFM,i:i-1+totNFM)=MM(1:totNFM,1:totNFM,n)
          n=n+1
       enddo
       A = B
       do i = 1,dim*ngauss
       write(100,'(10000f12.5)') (B(i,j), j=1,dim*ngauss)   
       enddo
end subroutine Matrix_A
!BL10
subroutine Flux_Guess(dim,ng,Nmat,nregion,ngauss,order,totNFM,NusigF,&
                      dcell,wt,p,phi_ni,phi_li)
       implicit none
       integer(kind=4), intent(in) :: totNFM,dim,ng,Nmat,nregion,ngauss,order
       real(kind=8), dimension(Nmat,ng), intent(in) :: NusigF
       real(kind=8), dimension(order,ngauss), intent(in) :: p
       real(kind=8), dimension(nregion), intent(in) :: dcell
       real(kind=8), dimension(ngauss), intent(in) :: wt
       real(kind=8), dimension(dim*ngauss), intent(out) :: phi_ni
       real(kind=8), dimension(dim,order), intent(out) :: phi_li
       integer(kind=4) :: i,j,k1,i1,n,l,m,k0,k2
       real(kind=8), dimension(ng*nregion) :: a10,a11
       real(kind=8), dimension(dim*ngauss) :: flux1
       real(kind=8), dimension(dim,order) :: flux2
       i = 1; flux1=0.0d0; flux2=0.0d0
       do k1=1,ng
          do K2=1,nregion
             a10(i) = NusigF(k2,k1)
             a11(i) = dcell(k2)
             i = i + 1
          enddo
       enddo 
       flux2(:,1)  = 1.0d0/dot_product(a10,a11)
       k0 = 1; k2 = 1; j = 1
       do k1 = 1,ng
          do n = 1,ngauss/2
             do i=k0,totNFM*k1
                l = 0
                do m = 1,order
                   flux1(j) = flux1(j) + 0.5d0*DBLE(2*l+1)*flux2(i,m)*p(m,n)
                   l = l+1
                enddo 
                   j = j+1
             enddo  
          enddo
          do n = ngauss/2+1,ngauss
             do i=totNFM*k1,k0,-1
                l = 0
                do m = 1,order
                   flux1(j) = flux1(j) + 0.5d0*DBLE(2*l+1)*flux2(i,m)*p(m,n)
                   l = l+1
                enddo
                   j = j+1 
             enddo  
          enddo
          k0 = k0 + totNFM
          k2 = k2 + ngauss
       enddo      
       k0 = 1; k2 = 1; j  = 1; i1  = 1
       do k1 = 1,ng
          do i = k0,totNFM*k1
             do n = 1,ngauss/2
                flux2(i,:) = flux2(i,:) + wt(n)*flux1(j)*p(:,n)
                j = j + totNFM
             enddo 
             i1 = i1 + 1
             j  = i1
          enddo
          j=totNFM*ngauss;  i1  = totNFM*ngauss
          do i = k0,totNFM*k1,k0
             do n = ngauss,ngauss/2+1,-1
                flux2(i,:) = flux2(i,:) + wt(n)*flux1(j)*p(:,n)
                j  = j - totNFM
             enddo 
             i1 = i1 - 1
             j  = i1
          enddo
          k0 = k0 + totNFM
          k2 = k2 + ngauss
       enddo
       phi_ni = flux1
       phi_li = flux2
end subroutine Flux_Guess
!BL11
subroutine Fission_Source(ng,dim,ngauss,order,totNFM,F,flux_li,p,k_eff,FQ_ni)
       implicit none
       integer(kind=4), intent(in) :: dim,ngauss,ng,totNFM,order
       real(kind=8), intent(in) :: k_eff 
       real(kind=8), dimension(dim ,dim), intent(in) :: F
       real(kind=8), dimension(order,ngauss), intent(in) :: p
       real(kind=8), dimension(dim,order), intent(in) :: flux_li
       real(kind=8), dimension(dim*ngauss), intent(out) :: FQ_ni
       real(kind=8), dimension(dim,order) :: Q_li,phi_li
       integer(kind=4) :: i,j,k0,k1,k2,l,n,m
       k0=1; k2=1; j = 1
       phi_li = flux_li
       FQ_ni(:) = 0.0
       if (order >= 2) then  
          do i = 2, order
             phi_li(:,i) = 0.0
          enddo
       endif

       Q_li(:,1) = matmul(F(:,:),phi_li(:,1))/k_eff

      do k1 = 1,ng
          do n = 1,ngauss/2
             do i=k0,totNFM*k1
                l = 0
                do m = 1,order
                   FQ_ni(j) = FQ_ni(j) + 0.5d0*DBLE(2*l+1)*Q_li(i,m)*p(m,n)
                   l = l+1
                enddo 
                   j = j+1
             enddo  
          enddo

          do n = ngauss/2+1,ngauss
             do i=totNFM*k1,k0,-1
                l = 0
                do m = 1,order
                   FQ_ni(j) = FQ_ni(j) + 0.5d0*DBLE(2*l+1)*Q_li(i,m)*p(m,n)
                   l = l+1
                enddo
                   j = j+1 
             enddo  
          enddo
          k0 = k0 + totNFM
          k2 = k2 + ngauss
       enddo

end subroutine Fission_Source
!BL12
subroutine Scattering_Source(ng,dim,ngauss,order,totNFM,D,U,L,phi_li,p,SQ_ni)
       implicit none
       integer(kind=4), intent(in) :: dim,ngauss,ng,totNFM,order
       real(kind=8), dimension(dim,dim,order), intent(in) :: L,U,D
       real(kind=8), dimension(order,ngauss), intent(in) :: p
       real(kind=8), dimension(dim,order), intent(in) :: phi_li
       real(kind=8), dimension(dim*ngauss), intent(out) :: SQ_ni
       real(kind=8), dimension(dim,order) :: Q_li
       integer(kind=4) :: i,j,k0,k1,k2,ii,ll,n,m
       SQ_ni(:) = 0.0
       do ii = 1,order
          Q_li(:,ii) = matmul(D(:,:,ii),phi_li(:,ii)) + &
                       matmul(L(:,:,ii),phi_li(:,ii)) + &
                       matmul(U(:,:,ii),phi_li(:,ii))
       enddo

       k2=1; k0=1; j=1
       do k1 = 1,ng
          do n = 1,ngauss/2
             do i=k0,totNFM*k1
                ll = 0
                do m = 1,order
                   SQ_ni(j) = SQ_ni(j) + 0.5d0*DBLE(2*ll+1)*Q_li(i,m)*p(m,n)
                   ll = ll+1
                enddo 
                   j = j+1
             enddo  
          enddo

          do n = ngauss/2+1,ngauss
             do i=totNFM*k1,k0,-1
                ll = 0
                do m = 1,order
                   SQ_ni(j) = SQ_ni(j) + 0.5d0*DBLE(2*ll+1)*Q_li(i,m)*p(m,n)
                   ll = ll+1
                enddo
                   j = j+1 
             enddo  
          enddo
          k0 = k0 + totNFM
          k2 = k2 + ngauss
       enddo
end subroutine Scattering_Source
!BL13
subroutine inverse(aa,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
       implicit none 
       INTEGER         ,INTENT(in)    :: n
       DOUBLE PRECISION,INTENT(in)  ::  aa(n,n)
       DOUBLE PRECISION,INTENT(out) ::  c(n,n)
       double precision L(n,n), U(n,n), a(n,n), b(n), d(n), x(n)
       double precision coeff
       integer i, j, k
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
       L=0.0;U=0.0;b=0.0;a=aa
! step 1: forward elimination
       do k=1, n-1
          do i=k+1,n
             coeff=a(i,k)/a(k,k)
             L(i,k) = coeff
             do j=k+1,n
                a(i,j) = a(i,j)-coeff*a(k,j)
             end do
          end do
       end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
       do i=1,n
         L(i,i) = 1.0
       end do
! U matrix is the upper triangular part of A
       do j=1,n
         do i=1,j
           U(i,j) = a(i,j)
         end do
       end do

! Step 3: compute columns of the inverse matrix C
       do k=1,n
         b(k)=1.0
         d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
         do i=2,n
           d(i)=b(i)
           do j=1,i-1
             d(i) = d(i) - L(i,j)*d(j)
           end do
         end do
! Step 3b: Solve Ux=d using the back substitution
         x(n)=d(n)/U(n,n)
         do i = n-1,1,-1
           x(i) = d(i)
           do j=n,i+1,-1
             x(i)=x(i)-U(i,j)*x(j)
           end do
           x(i) = x(i)/u(i,i)
         end do
! Step 3c: fill the solutions x(n) into column k of C
         do i=1,n
           c(i,k) = x(i)
         end do
         b(k)=0.0
       end do
       !write(*,'(16f7.3)') c
end subroutine inverse
!BL14
    subroutine GSL_SOR(ng,dim,Max_it,totNFM,ngauss,order,Nmat,iter,inter,eps,wt,mu,&
                       D,F,U,L,p,A,SigT,phi_ni,phi_li,Delta,k_eff,phi)
       implicit none
       ! global variables
       integer(kind=4), intent(in) :: ng,dim,Max_it,totNFM,ngauss,order,Nmat
       real(kind=8), dimension(dim*ngauss,dim*ngauss), intent(in) :: A
       real(kind=8), dimension(ngauss), intent(in) :: wt,mu
       real(kind=8), dimension(dim,dim,order), intent(in) :: D,U,L
       real(kind=8), dimension(dim,dim), intent(in) :: F
       real(kind=8), dimension(totNFM), intent(in) :: Delta
       real(kind=8), dimension(order,ngauss), intent(in) :: p
       real(kind=8), dimension(Nmat,ng), intent(in) :: SigT
       real(kind=8), intent(in) :: eps
       real(kind=8), dimension(dim*ngauss), intent(inout) :: phi_ni
       real(kind=8), dimension(dim,order), intent(inout) :: phi_li
       real(kind=8), dimension(dim), intent(out) :: phi
       real(kind=8), intent(out) :: k_eff
       integer(kind=4), intent(out) :: iter,inter
       ! variables locale
       real(kind=8), dimension(dim*ngauss) :: SQ_ni,FQ_ni,flux0,flux00,flux,vect0,vect1
       real(kind=8), dimension(dim*ngauss,dim*ngauss) :: ainv,dinv,INI,H,LL,ai
       real(kind=8), dimension(dim,order) :: phi_li0
       real(kind=8), dimension(ng) :: moy
       real(kind=8) :: omega= 1.,r=1.
       real(kind=8) ::  k_eff0=1.,err_k_eff,err_phi,Del,Sig,muu,eval1, eval2,err_flux,norme 
       integer(kind=4) :: i,j,k0,k1,i1,k2,n
       INI = 0.0;LL=0.;H=0.;iter=0;k_eff=1.;dinv=0.0
       err_k_eff = 1.2
       err_flux = 1.5
!      the positive flux condition
       Del = minval(Delta)
       Sig = minval(SigT)
       muu = minval(mu) 
       
       if (Del*Sig > 2.0d0*abs(muu)) then
          print*,'Failed the positive flux condition.'
          stop
          endif 
       !----------------------
       do i = 1,dim*ngauss
          dinv(i,i) = A(i,i)
          INI(i,i) = 1.0
       enddo
       !print*,'A'
       !write(*,'(8f6.2)') transpose(A) 
       !print*,'dinv'
       !write(*,'(8f6.2)') transpose(dinv)
       do i = 1,dim*ngauss
          do j = i,dim*ngauss
             if (i.ne.j) then
             LL(j,i) = -A(j,i)
             endif
          enddo
       enddo
       !print*,'L'
       !write(*,'(8f6.2)') transpose(LL)
       call inverse((dinv-r*LL),ainv,dim*ngauss)
       H = matmul(ainv,(1. - omega)*dinv+(omega-r)*LL) ! JOR: Jacobi overrelaxation method
       phi_li0  = phi_li
       flux0    = phi_ni  
       flux00   = phi_ni
       ! Outer Iteration
       do while ( err_k_eff >= eps .and. err_flux >= eps )
                iter = iter + 1
                call Fission_Source(ng,dim,ngauss,order,totNFM,F,phi_li,p,k_eff,FQ_ni)
                if (iter >= max_it) then
                print*,'Failed to converge.'
                stop
                endif  
                vect0 =  matmul(ainv,FQ_ni)
                flux = 0.0
                k_eff0 = k_eff
                err_phi = 1.
                !Inner Iteration
                do while ( err_phi >= 0.0001 )
                         call Scattering_Source(ng,dim,ngauss,order,totNFM,D,U,L,phi_li,p,SQ_ni) 
                         vect1  = matmul(ainv,SQ_ni) 
                         phi_ni = matmul(H,phi_ni) + omega*(vect0 + vect1)
                         ! The Legendre moments of the flux 
                         !-----------------------------------------------
                         phi_li = 0.0; k0 = 1; k2 = 1; j  = 1; i1  = 1
                         do k1 = 1,ng
                            do i = k0,totNFM*k1
                               do n = 1,ngauss/2
                                  phi_li(i,:) = phi_li(i,:) + wt(n)*phi_ni(j)*p(:,n)
                                  !print*,wt(n),phi_ni(j),p(:,n)
                                  j = j + totNFM
                               enddo 
                               i1 = i1 + 1
                               j  = i1
                            enddo
                            j=totNFM*ngauss;  i1  = totNFM*ngauss
                            do i = k0,totNFM*k1
                               do n = ngauss,ngauss/2+1,-1
                                  phi_li(i,:) = phi_li(i,:) + wt(n)*phi_ni(j)*p(:,n)
                                  j = j - totNFM
                               enddo 
                               i1 = i1 - 1
                               j  = i1
                            enddo
                         k0 = k0 + totNFM
                         k2 = k2 + ngauss
                         enddo
                         !write(*,'(16F6.2)') phi_li
                         !-----------------------------------------------
                         flux =  flux + phi_ni  
                         err_phi = maxval(abs(flux-flux0)/flux0)           
                         flux0 = flux 
                         inter = inter + 1  
                         norme = sqrt(sum(phi_ni*phi_ni))
                         phi_ni = phi_ni/norme
                enddo ! ending intern Iteration
             err_flux =  maxval(abs(phi_ni-flux00)/flux00)     
             eval1=sum(matmul(F, phi_li(:,1)))
             eval2=sum(matmul(F, phi_li0(:,1)))
             k_eff = k_eff*(eval1/eval2)
             err_k_eff =  abs(k_eff-k_eff0)/k_eff
             write(*,2000)iter,k_eff,err_k_eff
             phi_li0 = phi_li
             flux00   = phi_ni
       end do  ! ending ixtern Iteration 
       !-----------------------------------------------
       ! calculating the scalar flux
       k0 = 1; k2 = 1; j  = 1; i1  = 1
       do k1 = 1,ng
          do i = k0,totNFM*k1
              do n = 1,ngauss/2
                 phi(i) =  phi(i) + wt(n)*phi_ni(j)
                 j = j + totNFM
              enddo 
                 i1 = i1 + 1
                 j  = i1
          enddo
          j=totNFM*ngauss;  i1  = totNFM*ngauss
          do i = k0,totNFM*k1,k0
             do n = ngauss,ngauss/2+1,-1
                phi(i) =  phi(i) + wt(n)*phi_ni(j)
                j = j - totNFM
             enddo 
                i1 = i1 - 1
                j  = i1
          enddo
          k0 = k0 + totNFM
          k2 = k2 + ngauss
       enddo 
       !Normalized flux 
       !------------------------------
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
       2000 format(t3,"Iteration",i4,":",5x,"===>",5x,"keff =",F9.6,5x,"===>",5x,"res =",e10.3)
    end subroutine GSL_SOR
!BL13
subroutine conjgrad(ng,dim,Max_it,totNFM,ngauss,order,Nmat,iter,inter,eps,wt,mu,&
                         D,F,U,L,p,A,SigT,phi_ni,phi_li,Delta,k_eff,phi) 
 ! Solve Ax=b with OrthoMin algo
 !  use constants
  implicit none
       integer(kind=4), intent(in) :: ng,dim,Max_it,totNFM,ngauss,order,Nmat
       real(kind=8), dimension(dim*ngauss,dim*ngauss), intent(in) :: A
       real(kind=8), dimension(ngauss), intent(in) :: wt,mu
       real(kind=8), dimension(dim,dim,order), intent(in) :: D,U,L
       real(kind=8), dimension(dim,dim), intent(in) :: F
       real(kind=8), dimension(totNFM), intent(in) :: Delta
       real(kind=8), dimension(order,ngauss), intent(in) :: p
       real(kind=8), dimension(Nmat,ng), intent(in) :: SigT
       real(kind=8), intent(in) :: eps
       real(kind=8), dimension(dim*ngauss), intent(inout) :: phi_ni
       real(kind=8), dimension(dim,order), intent(inout) :: phi_li
       real(kind=8), dimension(dim), intent(out) :: phi
       real(kind=8), intent(out) :: k_eff
       integer(kind=4), intent(out) :: iter,inter

  real(kind=8), dimension(dim*ngauss) :: SQ_ni,FQ_ni,Ap
  real(kind=8), dimension(dim*ngauss) :: xi,vi,ri
  real(kind=8), dimension(totNFM) :: moy
  real(kind=8) :: di=1.,dip,wi,psi,d0,eval1, eval2=1.0,err_k_eff=1.0,norm
  integer(kind=4) :: i,j,k0,k1,i1,k2,n,k3
  ! ... initial state ...
  iter = 0;inter = 0
  k_eff = 1.0
  xi  = phi_ni  ! initial guess
  ! ... CG MINRES iteration ...
      do while ( err_k_eff >= eps )
     !.........
      call Fission_Source(ng,dim,ngauss,order,totNFM,F,phi_li,p,k_eff,FQ_ni)
                if (iter >= max_it) then
                print*,'Failed to converge.'
                stop
                endif   
     !.........
           do k3=1,10
           call Scattering_Source(ng,dim,ngauss,order,totNFM,D,U,L,phi_li,p,SQ_ni) 
           ri  = MATMUL(A,xi) - FQ_ni - SQ_ni 
           vi  = -ri 
           di  = DOT_PRODUCT(ri,ri)
           norm= sqrt(di)
           do while (norm>=1.e-2)
                    Ap = MATMUL(A,vi)
                    dip = DOT_PRODUCT(vi,Ap)
                    wi  = di/dip
                    ri  = ri + Ap*wi
                    xi  = xi + vi*wi

                    d0  = DOT_PRODUCT(ri,ri)
                    psi = d0/di
                    di  = d0
                    vi  = -ri+vi*psi
                    norm= sqrt(di)
                    print*,norm
           enddo
           !dip = DOT_PRODUCT(vi,MATMUL(A,vi))

                         !-----------------------------------------------
                         phi_li = 0.0; k0 = 1; k2 = 1; j  = 1; i1  = 1
                         do k1 = 1,ng
                            do i = k0,totNFM*k1
                               do n = 1,ngauss/2
                                  phi_li(i,:) = phi_li(i,:) + wt(n)*xi(j)*p(:,n)
                                  !print*,wt(n),phi_ni(j),p(:,n)
                                  j = j + totNFM
                               enddo 
                               i1 = i1 + 1
                               j  = i1
                            enddo
                            j=totNFM*ngauss;  i1  = totNFM*ngauss
                            do i = k0,totNFM*k1
                               do n = ngauss,ngauss/2+1,-1
                                  phi_li(i,:) = phi_li(i,:) + wt(n)*xi(j)*p(:,n)
                                  j = j - totNFM
                               enddo 
                               i1 = i1 - 1
                               j  = i1
                            enddo
                         k0 = k0 + totNFM
                         k2 = k2 + ngauss
                         enddo
                         !-----------------------------------------------
                         inter = inter + 1 
           enddo
            eval1 = sum(matmul(F, phi_li(:,1)))
             k_eff = k_eff*(eval1/eval2)
             err_k_eff =  abs(eval2-eval1)/eval2
             eval2  = eval1
      iter   = iter + 1
      write(*,2000)iter,k_eff,err_k_eff
      enddo
       !-----------------------------------------------
       ! calculating the scalar flux
       k0 = 1; k2 = 1; j  = 1; i1  = 1
       do k1 = 1,ng
          do i = k0,totNFM*k1
              do n = 1,ngauss/2
                 phi(i) =  phi(i) + wt(n)*xi(j)
                 j = j + totNFM
              enddo 
                 i1 = i1 + 1
                 j  = i1
          enddo
          j=totNFM*ngauss;  i1  = totNFM*ngauss
          do i = k0,totNFM*k1,k0
             do n = ngauss,ngauss/2+1,-1
                phi(i) =  phi(i) + wt(n)*xi(j)
                j = j - totNFM
             enddo 
                i1 = i1 - 1
                j  = i1
          enddo
          k0 = k0 + totNFM
          k2 = k2 + ngauss
       enddo 
! ... return final result ...
! ... Normalized flux ...
      moy = 0.0
      k1 = 1
      do k0 = 1,ng
         do i=k1,totNFM*k0
         moy(k0) = moy(k0) + xi(i)
         end do
         do i=k1,totNFM*k0
         phi(i) = xi(i)/(moy(k0)/totNFM)
         end do
         k1 = k1 + totNFM 
      enddo
      return
! ... Format ...
      2000 format(t3,"Iteration",i4,":",5x,"===>",5x,"keff =",F9.6,5x,&
                  "===>",5x,"res =",e10.3)
end subroutine conjgrad
!BL14
subroutine Eigenvalues(ng,dim,Max_it,totNFM,ngauss,order,Nmat,it,inter,eps,wt,mu,&
                                   D,F,U,L,p,C,SigT,phi_ni,flux_li,Delta,k_eff,phi)
       implicit none 
       integer(kind=4), intent(in) :: ng,dim,totNFM,ngauss,order,Nmat,Max_it
       real(kind=8), dimension(ngauss), intent(in) :: wt,mu
       real(kind=8), dimension(totNFM), intent(in) :: Delta
       real(kind=8), dimension(dim ,dim), intent(in) :: F
       real(kind=8), dimension(order,ngauss), intent(in) :: p
       real(kind=8), dimension(Nmat,ng), intent(in) :: SigT
       real(kind=8), dimension(dim,dim,order), intent(in) :: D,U,L
       real(kind=8), dimension(dim*ngauss,dim*ngauss), intent(in) :: C
       real(kind=8), intent(in) :: eps
       real(kind=8), dimension(dim*ngauss), intent(inout) :: phi_ni
       real(kind=8), dimension(dim,order), intent(inout) :: flux_li
       real(kind=8), dimension(dim), intent(out) :: phi
       real(kind=8), intent(out) :: k_eff
       integer(kind=4), intent(out) :: it,inter
!      variables locale
       real(kind=8), dimension(dim*ngauss) :: Q_ni,SQ_ni,flux0,flux00,flux,FQ_ni 
       real(kind=8), dimension(dim,order) :: flux_li0
       real(kind=8), dimension(ng) :: moy 
       real(kind=8) :: err_k_eff, err_phi, k_eff0, eval1, eval2,Del,Sig,muu,err_flux
       integer(kind=4) :: i,j,m,n,k0,k2,k1,i1
!      convergence parameters
       !print*,phi_ni
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
       
       if (Del*Sig > 2.0d0*abs(muu)) then
          print*,'Failed the positive flux condition.'
          stop
          endif 
       
       flux_li0 = flux_li
       flux0    = phi_ni  
       flux00   = phi_ni
       !print*,flux_li
       do while ( err_k_eff >= eps .and. err_phi >= eps )
                 flux0 =  phi_ni
                 call Fission_Source(ng,dim,ngauss,order,totNFM,F,flux_li,p,k_eff,FQ_ni)
                 if (it >= max_it) then
                 print*,'Failed to converge.'
                 stop
                 endif 
                 flux = 0.0
                 k_eff0 = k_eff
                 err_flux = 1.5
! Staring intern Iteration 
! ==============================================================================
             do while (  err_flux >= 0.0001 )
                  
                 call Scattering_Source(ng,dim,ngauss,order,totNFM,D,U,L,flux_li,p,SQ_ni)
                 Q_ni = SQ_ni + FQ_ni
                 !write(*,'(1f8.4)')  Q_ni
                 phi_ni  = matmul(c,Q_ni)
                 !do i=1,dim*ngauss
                 !write(*,'(1f8.4,1f8.4)')  phi_ni(i),Q_ni(i)
                 !enddo
     
                 flux_li = 0.0; k0 = 1; k2 = 1; j  = 1; i1  = 1
                 do k1 = 1,ng
                    do i = k0,totNFM*k1
                       do n = 1,ngauss/2
                           flux_li(i,:) = flux_li(i,:) + wt(n)*phi_ni(j)*p(:,n)
                           !print*,wt(n),phi_ni(j),p(:,n)
                           j = j + totNFM
                       enddo 
                           i1 = i1 + 1
                           j  = i1
                    enddo
                           j=totNFM*ngauss;  i1  = totNFM*ngauss
                    do i = k0,totNFM*k1
                       do n = ngauss,ngauss/2+1,-1
                           flux_li(i,:) = flux_li(i,:) + wt(n)*phi_ni(j)*p(:,n)
                           j = j - totNFM
                       enddo 
                           i1 = i1 - 1
                           j  = i1
                    enddo
                       k0 = k0 + totNFM
                       k2 = k2 + ngauss
                 enddo
                 flux =  flux + phi_ni  
                 err_flux = maxval(abs(flux-flux0)/flux0)
                 
                 !write(*,*) inter,maxval(flux), maxval(flux0) 
                 flux0 = flux 
                 inter = inter + 1 
             enddo
! ending intern Iteration           

             err_phi =  maxval(abs(phi_ni-flux00)/flux00)
             it = it + 1     
             eval1=sum(matmul(F, abs(flux_li(:,1))))
             eval2=sum(matmul(F, abs(flux_li0(:,1))))
             k_eff = k_eff*(eval1/eval2)
             err_k_eff =  abs(k_eff-k_eff0)/k_eff
             write(*,2000)it,k_eff,err_k_eff
             flux_li0 = flux_li
             flux00   = phi_ni
       end do      

! ==============================================================================
! calculating the scalar flux
                 k0 = 1; k2 = 1; j  = 1; i1  = 1
                 do k1 = 1,ng
                    do i = k0,totNFM*k1
                       do n = 1,ngauss/2
                           phi(i) =  phi(i) + wt(n)*phi_ni(j)
                           j = j + totNFM
                       enddo 
                           i1 = i1 + 1
                           j  = i1
                    enddo
                           j=totNFM*ngauss;  i1  = totNFM*ngauss
                    do i = k0,totNFM*k1,k0
                       do n = ngauss,ngauss/2+1,-1
                           phi(i) =  phi(i) + wt(n)*phi_ni(j)
                           j = j - totNFM
                       enddo 
                           i1 = i1 - 1
                           j  = i1
                    enddo
                       k0 = k0 + totNFM
                       k2 = k2 + ngauss
                 enddo   
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
       2000 format(t3,"Iteration",i4,":",5x,"===>",5x,"keff =",F9.6,5x,"===>",5x,"res =",e10.3)
end subroutine Eigenvalues
!BL16
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
        open (100,file='app/Output/OUTPUT_SN.TXT')
        write (100, FMT=* ) '********************************************************************************'
        write (100, FMT=* ) 'ERSN, UNIVERSITY ABDELMALEK ESSAADI FACULTY OF SCIENCES - TETOUAN, MOROCCO'
        write (100, FMT=* ) 'CODE  DEVELOPED  BY  MOHAMED  LAHDOUR,  PHD  STUDENT'
        write (100, FMT=* ) 'NTP-ERSN:        SN  DISCRETE  ORDINATES  METHOD'
        write (100, FMT=* ) 'VERSION NUMBER:  1.2'
        write (100, FMT=* ) 'VERSION DATE:    8  OTOBER  2018'
        write (100,3010) 'RAN ON:          ', start,'(H:M:S)'
        write (100, FMT=* ) '********************************************************************************'
        write (100, FMT=* ) '           ----------------------------------------------------------' 
        write (100, FMT=* ) '                     INPUT  PARAMETER - VALUES  FROM  INPUT'              
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'ENERGY GROUPS NUMBER:                    ',ng
        write (100, FMT=* ) 'REGIONS NUMBER:                          ',nregion
        write (100, FMT=* ) 'MATERIALS NUMBER:                        ',Nmat
        write (100,3040)    'SIZE OF EACH REGION [CM]:                ',dcell       
        write (100, FMT=* ) 'ANGULAR DISCRETIZATIONS:                 ',ngauss
        write (100, FMT=* ) 'ORDER LEGENDRE POLYNOMIAL:               ',order-1
        write (100, FMT=* ) 'TOTAL NUMBER OF FINE MESHES:             ',totNFM
        write (100,3050)    'CONVERGENCE CRITERION of KEFF AND FLUX:  ',eps
        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '                      CALCULATION  RUN-TIME  PARAMETERS  SN' 
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
!BL17
subroutine plot_flux(Delta,flux,totNFM,dim)
        implicit none
        integer(kind=4), intent(in) :: dim,totNFM
        real(kind=8), dimension(dim), intent(in) :: flux
        real(kind=8), dimension(totNFM), intent(in) :: Delta
        real(kind=8) :: som
        integer(kind=4) :: i,j
        open (10,file='app/Output/flux_sn.h')
        som = Delta(1)
        do i=1,totNFM
        write(10,*) som,(flux(i+j), j=0,dim-1,totNFM)  
        som = som + Delta(i)
        enddo
        close(10)
end subroutine plot_flux
!BL18
subroutine title1()       
       write(*,FMT='(/20(A/))') &
       '      ███╗   ██╗████████╗██████╗&
       &       ███████╗██████╗ ███████╗███╗   ██╗',&
       '      ████╗  ██║╚══██╔══╝██╔══██╗      ██&
       &╔════╝██╔══██╗██╔════╝████╗  ██║',&
       '      ██╔██╗ ██║   ██║   ██████╔╝█████╗█████╗&
       &  ██████╔╝███████╗██╔██╗ ██   ',&
       '      ██║╚██╗██║   ██║   ██╔═══╝ ╚════╝██╔══╝&
       &  ██╔══██╗╚════██║██║╚██╗██║',&
       '      ██║ ╚████║   ██║   ██║           ███████╗██║&
       &  ██║███████║██║ ╚████║',&
       '      ╚═╝  ╚═══╝   ╚═╝   ╚═╝           ╚══════╝&
       &╚═╝  ╚═╝╚══════╝╚═╝  ╚═══╝',&
         '______________________________________________________________________________'
       write(*,FMT=*) '                                                   Version Number: 1.2 '
       write(*,FMT=*) '     Copyright:      2015-2018 FS-Tetouan University Abdelmalk Essaadi '
       write ( *, FMT=* ) '. '
       write ( *, FMT=* ) '   FORTRAN90 version'  
       write ( *, FMT=* ) '   The Discrete Ordinates Method Sn'  
       write ( *, FMT=* ) '   Calculation of 1D Discrete Angle Domain'
       write ( *, FMT=* ) '   Slab 1D geometry'
       write ( *, FMT=* ) '. '
end subroutine title1
!BL19
subroutine title2()
       write ( *, FMT=* )' ************************************************************************'
       write ( *, FMT=* )'                               Finished'                             
       write ( *, FMT=* )' ************************************************************************'  
end subroutine title2
!BL20
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

