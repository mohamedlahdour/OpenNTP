!>BL1
subroutine gauleg(x1,x2,x,w,ngauss)
       !> EPS is the relative precision.
       !> Given the lower and upper limits of integration x1 and x2 ,
       !> and given n , this routine returns
       !> arrays x(1:n) and w(1:n) of length n ,
       !> containing the abscissas and weights of the Gauss-
       !> Legendre n-point quadrature formula.
       implicit none
       integer(kind=4), intent(in) :: ngauss
       real(kind=4), intent(in) :: x1,x2
       real(kind=4), dimension(ngauss), intent(out) :: x,w
       real(kind=4), parameter :: EPSS = 3.0e-7
       integer(kind=4) :: i,j,m
       real(kind=4) :: p1,p2,p3,pp,xl,xm,z,z1
       m=(ngauss+1)/2
       xm=0.50*(x2+x1)
       xl=0.50*(x2-x1)
       do  i=1,m         !Loop over the desired roots.
       z=cos(3.141592654*(i-.25)/(ngauss+.50))
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
!>BL2
subroutine leg_poly(ngauss,order,mu,p) 
       !> p  : Legendre Polynomials
       !> mu : direction cosine
       !> ngauss-point quadrature formula
       !> order : The l-order Legendre polynomial
       implicit none
       integer(kind=4), intent(in) :: ngauss,order
       real(kind=4), dimension(ngauss), intent(in) :: mu
       real(kind=4), dimension(order,ngauss), intent(out) :: p
       real(kind=4) :: a1,a2
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
!>BL3
    subroutine fmm_id(assembly,core,fmmid,nfmesh,width,RegMat,npx,npc,nx,nxx,na,totNFM,Delta)
       !>  This subroutine fmm_id(i) allowing to identify which type of material in the cell i.
       !>  npx     : Total number of Coarse mesh along x axis for each pincell.
       !>  npc     : Total number of PinCell
       !>  na      : Total number of assembly
       !>  nxx     : Total number of Coarse mesh along x axis for each assembly.
       !>  nx      : Total number of Coarse mesh along x axis.
       !>  nfmesh  : Numbre of finite mesh in each pin cell along npx axis.
       !>  width   : Number of fine mesh per x coarse mesh for each pincell.
       !>  core    : reactor core
       implicit none
       integer(kind=4), intent(in) :: npx,npc,nx,nxx,na,totNFM
       integer(kind=4), dimension(npc,npx), intent(in) :: RegMat, nfmesh
       real(kind=4), dimension(npc,npx), intent(in) :: width
       integer(kind=4), dimension(nx), intent(in) :: core
       integer(kind=4), dimension(na,nxx), intent(in) :: assembly
       integer(kind=4), dimension(totNFM), intent(out) :: fmmid
       real(kind=4), dimension(totNFM), intent(out) :: Delta
       ! Variables Locales
       integer(kind=4) :: i,n1,n2,n3
       i = 1
       ! -- fine mesh material id 
       do n1 = 1,nx
          do n2 = 1,nxx
             do n3 = 1,npx    
                fmmid(i:i+nfmesh(assembly(core(n1),n2),n3)-1) = RegMat(assembly(core(n1),n2),n3)
                Delta(i:i+nfmesh(assembly(core(n1),n2),n3)-1) = width(assembly(core(n1),n2),n3)/&
                                                                nfmesh(assembly(core(n1),n2),n3)
                i=i+nfmesh(assembly(core(n1),n2),n3)
             enddo
          enddo
       enddo
    end subroutine fmm_id
!>BL5
subroutine Matrix_D(SigS,D,fmmid,ng,Nmat,order,totNFM,dim)
       !> Diagonal D matrix contains scattering cross section in each cell
       !> SigS   : Scattering cross section 
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,totNFM,dim,order
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=4), dimension(Nmat,order,ng,ng), intent(in) :: SigS
       real(kind=4), dimension(dim,dim,order), intent(out) :: D
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
!>BL6
subroutine Matrix_L(SigS,L,fmmid,ng,Nmat,order,totNFM,dim)
       !> Lower matrix of the scattering cross section
       !> SigS : Scatternig cross section
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,totNFM,dim,order
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=4), dimension(Nmat,order,ng,ng), intent(in) :: SigS
       real(kind=4), dimension(dim ,dim,order), intent(out) :: L
       integer(kind=4) :: i,k0,k1,k2,k3,k4
                L(:,:,:) = 0.0
            do  k4  = 1,order
                k2  = 0
            do  k1  = 1,ng
                k0  = 1
                k3  = 0
                    do while (k0<k1) 
                             do  i  = 1,totNFM
                             L(i+k3,i+(k1-1)*totNFM,k4) = SigS(fmmid(i),k4,k1,k0) ! 1 --> 2 SigS(1,1,2)
                             enddo
                             k0 = k0 + 1
                             k3 = k3 + totNFM
                    enddo
            enddo
            enddo
end subroutine Matrix_L
!>BL7
subroutine Matrix_U(SigS,U,fmmid,ng,Nmat,order,totNFM,dim)
       !> Upper matrix of the scattering cross section
       !> SigS : Scatternig cross section
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,totNFM,dim,order
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=4), dimension(Nmat,order,ng,ng), intent(in) :: SigS
       real(kind=4), dimension(dim,dim,order), intent(out) :: U
       integer(kind=4) :: i,k0,k1,k2,k3,k4
                U(:,:,:) = 0.0
            do  k4  = 1,order
                k2  = 0
            do  k1  = 1,ng
                k0  = 1
                k3  = 0
                    do while (k0<k1) 
                             do  i  = 1,totNFM
                             U(i+(k1-1)*totNFM,i+k3,k4) = SigS(fmmid(i),k4,k0,k1) ! 1 --> 2 SigS(1,1,2)
                             enddo
                             k0 = k0 + 1
                             k3 = k3 + totNFM
                    enddo
            enddo
            enddo
end subroutine Matrix_U
!>BL8
subroutine Matrix_F(NusigF,Chi,F,fmmid,ng,Nmat,totNFM,dim)
       !> Matrix of the NuFission cross section
       !> NuSigF : NuFission cross section
       !> CHI    : Density function for neutrons
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,totNFM,dim
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=4), dimension(Nmat,ng), intent(in) :: NusigF,Chi
       real(kind=4), dimension(dim ,dim), intent(out) :: F
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
                             F(i+k3,i+(k1-1)*totNFM) =  Chi(fmmid(i),k0)*NusigF(fmmid(i),k1)
                             F(i+(k1-1)*totNFM,i+k3) =  Chi(fmmid(i),k1)*NusigF(fmmid(i),k0)
                             enddo
                             k0 = k0 + 1
                             k3 = k3 + totNFM
                    enddo
            enddo
end subroutine Matrix_F
!BL9
subroutine Matrix_AB(ng,Nmat,dim,totNFM,ngauss,mu,fmmid,SigT,Delta,A,B)
       !> construct matrix A and B
       !> A and B are non-symstric matrices resulting from the discretitation 
       !> of the transport equation
       !> generally A and B presents a block structure by energy group
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,dim,totNFM,ngauss
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=4), dimension(totNFM), intent(in) :: Delta
       real(kind=4), dimension(ngauss), intent(in) :: mu
       real(kind=4), dimension(Nmat,ng), intent(in) :: SigT
       real(kind=4), dimension(dim,ngauss*ng), intent(out) :: A,B
       integer(kind=4) :: i,n,k0,k1,k2,j
       real(kind=4) :: denom
       k0=1;k2=1
     
       do k1=1,ng
                 j = totNFM
             do  i = totNFM*k1,k0,-1 ! right-to-left 
                 do  n = 1,(ngauss/2)
                     denom  = 2*mu(n)-SigT(fmmid(j),k1)*Delta(j)
                     A(i,n + k2 -1) = (2*mu(n)+SigT(fmmid(j),k1)*Delta(j))/denom
                     B(i,n + k2 -1) = 2*Delta(j)/denom
                 enddo
                 j = j-1
             enddo

                     
                 j = 1
             do  i = k0,totNFM*k1    ! left-to-rights
                 do  n = (ngauss/2)+1,ngauss
                     denom  = 2*mu(n)+SigT(fmmid(j),k1)*Delta(j)
                     A(i,n + k2 -1) = (2*mu(n)-SigT(fmmid(j),k1)*Delta(j))/denom
                     B(i,n + k2 -1) = 2*Delta(j)/denom
                 enddo
                     j = j + 1
             enddo

             k0 = k0 + totNFM
             k2 = k2 + ngauss
       enddo
end subroutine Matrix_AB
!>BL10
subroutine flux_guess(dim,ng,Nmat,ngauss,order,totNFM,fmmid,NusigF,Delta,p,flux_ni,flux_li)
       !> This function allows to calculate the guess flux
       !> flux_ni   : The guess angular flux
       !> flux_li   : The Legendre moments of the guess flux 
       !> Delta     : The volume of each cell
       !> NusigF    : NuFission cross section
       !> p         : Legendre Polynomials
       implicit none
       integer(kind=4), intent(in) :: dim,ng,Nmat,ngauss,order,totNFM
       real(kind=4), dimension(Nmat,ng), intent(in) :: NusigF
       real(kind=4), dimension(order,ngauss), intent(in) :: p
       real(kind=4), dimension(totNFM), intent(in) :: Delta
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=4), dimension(dim,ngauss*ng), intent(out) :: flux_ni
       real(kind=4), dimension(dim,order), intent(out) :: flux_li
       integer(kind=4) :: i,k1,n,ll,m,k0,k2
       real(kind=4), dimension(ng*totNFM) :: a10,a11
       i = 1
       do k1=1,ng
          do k2=1,totNFM
             a10(i) = NusigF(fmmid(k2),k1)
             a11(i) = Delta(k2)
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
!>BL11
subroutine Fission_Source(ng,dim,ngauss,order,totNFM,F,flux_li,p,k_eff,Q_li,FQ_ni)
       !> Calculate the fission source components FQ_ni
       !> dim = totNFM x ng
       !> totNFM  : The total number of meshes
       !> ng      : The total number of energy group
       !> ngauss  : Number of angular discretizations
       !> flux_li : The Legendre moments of the flux 
       !> Q_li    : The Legendre moments of the fission source 
       !> p       : Legendre Polynomials
       implicit none
       integer(kind=4), intent(in) :: dim,ngauss,ng,totNFM,order
       real(kind=4), intent(in) :: k_eff 
       real(kind=4), dimension(dim ,dim), intent(in) :: F
       real(kind=4), dimension(order,ngauss), intent(in) :: p
       real(kind=4), dimension(dim,order), intent(in) :: flux_li
       real(kind=4), dimension(dim,ngauss*ng), intent(out) :: FQ_ni
       real(kind=4), dimension(dim,order), intent(out) :: Q_li
       real(kind=4), dimension(dim,order) :: phi_li
       integer(kind=4) :: i,k0,k1,k2,ll,n,m
       k0=1;k2=1
       phi_li = flux_li
       FQ_ni(:,:) = 0.0
       if (order >= 2) then  
          do i = 2, order
             phi_li(:,i) = 0.0
          enddo
       endif

       Q_li(:,:) = matmul(F(:,:),phi_li(:,:))/k_eff

       do k1 = 1,ng
          do i =k0,totNFM*k1
              do n = 1,ngauss
                ll = 0
                do m = 1,order
                    FQ_ni(i,n+k2-1) = FQ_ni(i,n+k2-1) + 0.5*(2.*float(ll)+1.)*Q_li(i,m)*p(m,n)
                   ll = ll+1
                enddo
              enddo   
          enddo
             k2 = k2 + ngauss
             k0 = k0 + totNFM
       enddo
end subroutine Fission_Source
!>BL12
subroutine Scattering_Source(ng,dim,ngauss,order,totNFM,D,U,L,flux_li,p,SQ_ni)
       !> Calculate the scattering source components SQ_ni
       !> dim = totNFM x ng
       !> totNFM  : The total number of meshes
       !> ng      : The total number of energy group
       !> ngauss  : Number of angular discretizations
       !> flux_li : The Legendre moments of the flux 
       !> Q_li    : The Legendre moments of the fission source 
       !> p       : Legendre Polynomials
       implicit none
       integer(kind=4), intent(in) :: dim,ngauss,ng,totNFM,order
       real(kind=4), dimension(dim,dim,order), intent(in) :: L,U,D
       real(kind=4), dimension(order,ngauss), intent(in) :: p
       real(kind=4), dimension(dim,order), intent(in) :: flux_li
       real(kind=4), dimension(dim,ngauss*ng), intent(out) :: SQ_ni
       real(kind=4), dimension(dim,order) :: Q_li
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
!>BL13
subroutine Total_Source(ng,dim,ngauss,FQ_ni,SQ_ni,Q_ni)
       !> Calculate the Total source components Q_ni
       !> dim = totNFM x ng
       !> totNFM  : The total number of meshes
       !> ng      : The total number of energy group
       !> ngauss  : Number of angular discretizations
       !> FQ_ni   : The fission source components
       !> SQ_ni   : The scattering source components
       implicit none
       integer(kind=4), intent(in) :: dim,ngauss,ng
       real(kind=4), dimension(dim,ngauss*ng), intent(in) :: FQ_ni,SQ_ni
       real(kind=4), dimension(dim,ngauss*ng), intent(out) :: Q_ni    
       Q_ni =  SQ_ni + FQ_ni
end subroutine Total_Source
!BL14
subroutine current_f(BC,A,B,Q_ni,curr,ng,ngauss,dim,totNFM)
       !> Sweep over a 1D mesh.
       !> BC    : Boundary conditions
       !> Q_ni  : Total source components
       !> curr  : The flux components in interfaces
       implicit none
       CHARACTER(50), intent(in) :: BC
       integer(kind=4), intent(in) :: ng,ngauss,dim,totNFM
       real(kind=4), dimension(dim,ngauss*ng), intent(in) :: A,B,Q_ni
       real(kind=4), dimension((totNFM+1)*ng,ngauss*ng), intent(out) :: curr
       real(kind=4), dimension((totNFM+1)*ng,ngauss*ng) ::curra,currb
       integer(kind=4) :: i,j,k1,k0,k2,m,n
       curr = 0.0
       k0=1;k2=1
       !========================\\reflective//==============================
        if  ( BC == 'Reflective' ) then
            do k1 = 1,ng
                        m = k1*totNFM
                do  i = k1*totNFM+k1-1,k0,-1 ! right-to-left 
                    do  j=1,(ngauss/2)
                        curra(k1*totNFM+k1,j + k2 -1) = 0.0
                        currb(k1*totNFM+k1,j + k2 -1) = 1.0
                        curra(i,j+k2-1) = A(m,j+k2-1)*curra(i+1,j+k2-1) - B(m,j+k2-1)*Q_ni(m,j+k2-1)
                        currb(i,j+k2-1) = A(m,j+k2-1)*currb(i+1,j+k2-1) - B(m,j+k2-1)*Q_ni(m,j+k2-1)
                        !   Calculate the current in the cell
                        curr(i+1,j+k2-1)= curra(i,j+k2-1)/ (1.0 + (curra(i,j+k2-1)-currb(i,j+k2-1)))
                    enddo
                        m = m -1
                enddo  

                        n = m + 1
                do  i = k0,k1*totNFM+k1-1    ! left-to-right
                    do  j=ngauss,(ngauss/2)+1,-1
                        curra(k0,j+k2-1) = 0.0
                        currb(k0,j+k2-1) = 1.0
                        curra(i+1,j+k2-1) = A(n,j+k2-1)*curra(i,j+k2-1) + B(n,j+k2-1)*Q_ni(n,j+k2-1)
                        currb(i+1,j+k2-1) = A(n,j+k2-1)*currb(i,j+k2-1) + B(n,j+k2-1)*Q_ni(n,j+k2-1)
                        !   Calculate the current in the cell
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
                        curr(i,j+k2-1) = A(m,j+k2-1)*curr(i+1,j+k2-1) - B(m,j+k2-1)*Q_ni(m,j+k2-1)
                    enddo
                        m = m -1
                enddo  
                        n = m + 1
                do  i = k0,k1*totNFM+k1-1    ! left-to-right
                    do  j=ngauss,(ngauss/2)+1,-1
                        curr(i+1,j+k2-1) = A(n,j+k2-1)*curr(i,j+k2-1) + B(n,j+k2-1)*Q_ni(n,j+k2-1) 
                    enddo
                        n = n + 1
                enddo
                        k0 = k0 + (totNFM + 1)
                        k2 = k2 + ngauss
           enddo
       !========================\\vacuum//==============================

        elseif  ( BC == 'Vacuum' ) then 
            do k1 = 1,ng
                        m = k1*totNFM
                do  i = k1*totNFM+k1-1,k0,-1 ! right-to-left 
                    do  j=1,(ngauss/2)
                        curr(k1*totNFM+k1,j + k2 -1) = 0.0
                        curr(i,j+k2-1) = A(m,j+k2-1)*curr(i+1,j+k2-1) - B(m,j+k2-1)*Q_ni(m,j+k2-1)
                    enddo
                        m = m -1
                enddo
                    
                        n = m + 1
                do  i = k0,k1*totNFM+k1-1    ! left-to-right
                    do  j=ngauss,(ngauss/2)+1,-1
                        curr(k0,j + k2 -1) = 0.0
                        curr(i+1,j+k2-1) = A(n,j+k2-1)*curr(i,j+k2-1) + B(n,j+k2-1)*Q_ni(n,j+k2-1) 
                    enddo
                        n = n + 1
                enddo
                    k0 = k0 + (totNFM + 1)
                    k2 = k2 + ngauss
            enddo
       !========================\\vacuum_reflective//====================
        elseif  ( BC == 'Vacuum Reflective' ) then
             do k1 = 1,ng
                    n = (k1-1)*totNFM + 1
                do  i = k0,k1*totNFM+k1-1    ! left-to-right
                    do  j=(ngauss/2+1),ngauss
                        curr(k0,j + k2 -1) = 0.0
                        curr(i+1,j+k2-1) = A(n,j+k2-1)*curr(i,j+k2-1) + B(n,j+k2-1)*Q_ni(n,j+k2-1) 
                    enddo
                    n = n + 1
                enddo
   
                    m = n - 1
                do  i = k1*totNFM+k1-1,k0,-1 ! right-to-left
                    do  j=1,(ngauss/2)
                        curr(k1*totNFM+k1,j+k2 -1) = curr(k1*totNFM+k1,ngauss-j+k2)
                        curr(i,j+k2-1) = A(m,j+k2-1)*curr(i+1,j+k2-1) - B(m,j+k2-1)*Q_ni(m,j+k2-1)
                    enddo
                    m = m - 1
                enddo 
                k0 = k0 + (totNFM + 1)
                k2 = k2 + ngauss
             enddo
       !========================\\reflective_vacuum//====================
        elseif  ( BC == 'Reflective Vacuum' ) then
             do k1 = 1,ng
                    m = k1*totNFM
                do  i = k1*totNFM+k1-1,k0,-1 ! right-to-left 
                    do  j=(ngauss/2),1,-1
                        curr(k1*totNFM+k1,j + k2 -1) = 0.0
                        curr(i,j+k2-1) = A(m,j+k2-1)*curr(i+1,j+k2-1) - B(m,j+k2-1)*Q_ni(m,j+k2-1)
                    enddo
                    m = m -1
                enddo
                    n = m + 1
                do  i = k0,k1*totNFM+k1-1    ! left-to-right
                    do  j=(ngauss/2+1),ngauss
                        curr(k0,j + k2 -1) = curr(k0,ngauss-j+k2) 
                        curr(i+1,j+k2-1) = A(n,j+k2-1)*curr(i,j+k2-1) + B(n,j+k2-1)*Q_ni(n,j+k2-1)
                    enddo
                    n = n + 1
                enddo
                k0 = k0 + (totNFM + 1)
                k2 = k2 + ngauss
             enddo
        endif
end subroutine current_f
!>BL15
subroutine Outer_Iteration(ng,dim,Max_it,totNFM,ngauss,order,Nmat,it,inter,scheme,eps,wt,mu,&
                          D,F,U,L,A,B,p,BC,SigT,Delta,k_eff,phi,NusigF,sigF,Chi,fmmid)
       !> Iterative method to calculate eigenvalues and vectors
       !> phi       : Scalar flux
       !> Max_it    : Maximum number of iterations
       !> k_eff     : multiplication factor
       !> eps       : Criterion of convergence
       !> the weights wt and base points mu are chose of the classical Gauss-Legendre quadrature. 
       implicit none 
       integer(kind=4), intent(in) :: ng,dim,totNFM,ngauss,order,Nmat,Max_it
       real(kind=4), dimension(ngauss), intent(in) :: wt,mu
       real(kind=4), dimension(totNFM), intent(in) :: Delta
       real(kind=4), dimension(dim,ngauss*ng), intent(in) :: A,B
       real(kind=4), dimension(dim ,dim), intent(in) :: F
       real(kind=4), dimension(order,ngauss), intent(in) :: p
       real(kind=4), dimension(Nmat,ng), intent(in) :: SigT
       real(kind=4), dimension(dim,dim,order), intent(in) :: D,U,L
       real(kind=4), intent(in) :: eps
       real(kind=4), dimension(Nmat,ng), intent(in) :: NusigF,sigF,Chi
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=4), dimension(dim), intent(out) :: phi
       real(kind=4), intent(out) :: k_eff
       integer(kind=4), intent(out) :: it,inter
!      local variables
       real(kind=4), dimension(dim,ngauss*ng) :: flux_ni
       real(kind=4), dimension(dim,order) :: flux_li,Q_li,flux_li0,flux_li1
       real(kind=4), dimension(dim,ngauss*ng) :: Q_ni,SQ_ni,FQ_ni 
       real(kind=4), dimension((totNFM+1)*ng,ngauss*ng) :: curr
       real(kind=4), dimension(ng) :: moy 
       real(kind=4) :: err_k_eff,err_phi,k_eff0,Del,Sig,muu,err_flux,alpha
       real(kind=4) :: dsnew,dsold,err1,err2
       integer(kind=4) :: i,j,m,n,k0,k2,k1
       CHARACTER(50), intent(in) :: BC, scheme

       call flux_guess(dim,ng,Nmat,ngauss,order,totNFM,fmmid,NusigF,Delta,p,flux_ni,flux_li)
       k_eff     = 1.0; k_eff0    = 1.0
       call Fission_Source(ng,dim,ngauss,order,totNFM,F,flux_li,p,k_eff,Q_li,FQ_ni)
       dsold = sum(Q_li)
!      convergence parameters
       err_k_eff = 1.0
       err_phi = 1.0
       err_flux = 1.5
       it = 0 
!      the positive flux condition
       Del = minval(Delta)
       Sig = minval(SigT)
       muu = minval(mu) 
       
       if (Del*Sig > 2*abs(muu)) then
          print*,'Failed the positive flux condition.'
          stop
          endif        
       do while ( err_k_eff >= eps .and. err_phi >= eps )
                 if (it >= max_it) then
                 print*,'Failed to converge.'
                 stop
                 endif 
                 flux_li1 = flux_li
                 err_flux = 1.0
                 inter=0
! Staring intern Iteration 
! ==============================================================================
             do  while (  err_flux >= eps*10)
                 flux_li0 = flux_li
                 call Scattering_Source(ng,dim,ngauss,order,totNFM,D,U,L,flux_li,p,SQ_ni)  
                 call Total_Source(ng,dim,ngauss,FQ_ni,SQ_ni,Q_ni)
                 call current_f(BC,A,B,Q_ni,curr,ng,ngauss,dim,totNFM)
                 k0 = 1
                 k2 = 1
            do k1 = 1,ng
                        m = k1*totNFM+k1-1
                do  i = k1*totNFM,k0,-1 ! right-to-left 
                    do  j=(ngauss/2),1,-1
                        if (scheme == 'Diamond Difference') then
                           alpha = 0
                        else
                           alpha = abs(mu(j))/mu(j)
                        endif
                        flux_ni(i,j+k2-1) = 0.5*((1+alpha)*curr(m+1,j+k2-1) + (1-alpha)*curr(m,j+k2-1))
                    enddo
                        m = m - 1
                enddo        
                        n = m + 1
                do  i = k0,k1*totNFM    ! left-to-right
                    do  j=ngauss,(ngauss/2)+1,-1
                        if (scheme == 'Diamond Difference') then
                           alpha = 0
                        else
                           alpha = abs(mu(j))/mu(j)
                        endif
                        flux_ni(i,j+k2-1) = 0.5*((1+alpha)*curr(n+1,j+k2-1) + (1-alpha)*curr(n,j+k2-1)) 
                    enddo
                        n = n + 1
                enddo
                        k0 = k0 + totNFM 
                        k2 = k2 + ngauss
            enddo

                 flux_li = 0.0
                 k0 = 1
                 k2 = 1
                 do k1 = 1,ng
                    do i  = k1*totNFM,k0,-1
                       do n = 1,(ngauss/2)
                          flux_li(i,:) = flux_li(i,:) + wt(n)*flux_ni(i,n+k2-1)*p(:,n)
                       enddo
                    enddo 

                    do i = k0,totNFM*k1
                       do n = (ngauss/2)+1,ngauss
                           flux_li(i,:) = flux_li(i,:) + wt(n)*flux_ni(i,n+k2-1)*p(:,n)
                       enddo
                    enddo
                       k0 = k0 + totNFM
                       k2 = k2 + ngauss
                 enddo
                !Condition on the scalar flux
                err1   = maxval(abs(flux_li0))
                err2   = maxval(abs(flux_li-flux_li0))
                err_flux = err2/err1 
                !print*, err_flux 
                inter = inter + 1 
                if (inter>1000) exit 
             enddo
! ending intern Iteration 
! ==============================================================================
             ! Condition on the scalar flux
             err1   = maxval(abs(flux_li1))
             err2   = maxval(abs(flux_li-flux_li1))
             err_phi =  err2/err1
             call Fission_Source(ng,dim,ngauss,order,totNFM,F,flux_li,p,k_eff,Q_li,FQ_ni)
             !Normalised Source
             call NormalizeFlux(dim,totNFM,Nmat,ng,sigF,fmmid,delta,flux_li(:,1)) 
             dsnew = sum(Q_li)
             k_eff =  k_eff*dsnew/dsold
             err_k_eff =  abs(k_eff-k_eff0)/k_eff0
             it = it + 1
             write(*,2000)it,k_eff,err_k_eff     
             dsold = dsnew
             k_eff0   = k_eff
             if (it>500) exit 
       end do  ! ending extern Iteration  
             phi = flux_li(:,1)
       2000 format("Iteration",i4,":",1x,"===>",1x,"keff =",F9.6,1x,"===>",1x,"res =",e10.3)
end subroutine Outer_Iteration
!>BL16
subroutine Output(start,tm,k_eff,SigT,NusigF,SigS,Chi,mu,wt,xcm,phi,eps,TNFM_x,dim,ng,&
                  Nmat,order,npx,N,it1,it2,npc,na,nx,nxx,SFPC,PF)
        implicit none
        ! Global variables
        integer(kind=4), intent(in) :: dim,ng,TNFM_x,Nmat,order,npx,N,it1,it2,npc,na,nx,nxx
        real(kind=4), dimension(Nmat,ng), intent(in) :: SigT,NusigF,Chi
        real(kind=4), dimension(Nmat,order,ng,ng), intent(in) :: SigS
        real(kind=4), dimension(dim), intent(in) :: phi
        real(kind=4), dimension(npc,npx), intent(in)  :: xcm
        real(kind=4), dimension(nx*nxx,ng), intent(in) :: SFPC
        real(kind=4), dimension(nx*nxx), intent(in) :: PF
        real(kind=4), dimension(N), intent(in) :: mu,wt
        CHARACTER(50), intent(in) :: start,tm
        real(kind=4), intent(in) :: eps,k_eff
        ! Local variables
        real(kind=4), dimension(TNFM_x,ng) :: flux
        integer(kind=4) :: i,j,k,m
        m=1
        do i=1,ng
           do j=1,TNFM_x
              flux(j,i)=phi(m)
              m=m+1
           enddo
        enddo
        open (100,file='app/Output/OUTPUT_SN.TXT')
        write (100, FMT=* ) '********************************************************************************'
        write (100, FMT=* ) 'ERSN, UNIVERSITY ABDELMALEK ESSAADI FACULTY OF SCIENCES - TETOUAN, MOROCCO'
        write (100, FMT=* ) 'CODE  DEVELOPED  BY  MOHAMED  LAHDOUR,  PHD  STUDENT'
        write (100, FMT=* ) 'OpenNTP:         SN  DISCRETE  ORDINATES  METHOD'
        write (100, FMT=* ) 'DIMENSION:       ONE DIMENSION (1D) '
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
        write (100, FMT=* ) 'MATERIALS NUMBER:                        ',Nmat
        do i=1,npc
        write (100,3040)    'SIZE OF EACH X REGION PIN CELL',i,':     ',xcm(i,:)   
        enddo  
        write (100, FMT=* ) 'DISCRETIZATIONS ANGULAR:                 ',N
        write (100, FMT=* ) 'ORDER LEGENDRE POLYNOMIAL:               ',order-1
        write (100, FMT=* ) 'TOTAL NUMBER OF X FINE MESHES:           ',TNFM_x
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
        do i=1,N
        write(100,3060) i,mu(i),wt(i)
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
            write(100,3081) '|',j, sum(SigT(i,j)*flux(:,j)), sum((SigT(i,j)-SigS(i,1,j,j))*flux(:,j)),&
                               sum(NusigF(i,j)*flux(:,j)),sum(SigS(i,1,j,j)*flux(:,j)),'    |'
            enddo
        enddo
        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '                              NORMALIZED PIN POWER                   '   
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write(100,'(1x,2000A13)') ('-------------', i=1,nx*nxx+1)
        write(100,'(1x,A1,i14,2000i13)') '|',(i, i=1,nx*nxx)
        write(100,'(1x,2000A13)') ('-------------', i=1,nx*nxx+1)
        j=1
        write(100,2000) '|',j,(PF(i), i=1,nx*nxx)
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
        write(100,'(1x,2000A13)') ('-------------', i=1,nx*nxx+1)
        write(100,'(1x,A1,i14,2000i13)') '|',(i, i=1,nx*nxx)
        write(100,'(1x,2000A13)') ('-------------', i=1,nx*nxx+1)
        j=1
        write(100,2000) '|',j,(SFPC(i,k), i=1,nx*nxx) 
        enddo
        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '                       NORMALIZED SCALAR  FLUX  SOLUTION             ' 
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'FLUXES  PER  MESH  PER  ENERGY  GROUP:' 
        write (100, FMT=* ) '' 
        write (100,6000)'       M E S H ', ('     G R O U P',i,i=1,ng)
        write (100, FMT=* ) ''
        do i = 1,TNFM_x
        write(100,5000) i,(flux(i,j),j=1,ng)
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
        3060 format(1x,1p,i11,5x,e16.5,e16.5)
        3070 format(1x,A18,i4)
        3080 format(1x,1p,A1,i11,5x,e16.5,e16.5,e16.5,e16.5,e16.5,A12)
        3081 format(1x,1p,A1,i11,5x,e16.5,e16.5,e16.5,e16.5,A7)
        3090 format(1x,A26,6x,f8.6)
        4000 format(1x,A26,4x,A10,A10)
        5000 format(1x,1p,i11,5x,200e16.5)
        6000 format(1x,A14,2x,300(A14,i2)) 
        close(100)
end subroutine Output
!>BL17
subroutine Plot_flux(dim,totNFM,Nmat,ng,nx,nxx,napc,Delta,phi,SFPC,SF)
        !> Plot scalar flux
        implicit none
        integer(kind=4), intent(in) :: dim,totNFM,Nmat,ng,nx,nxx,napc
        real(kind=4), dimension(totNFM), intent(in) :: delta
        real(kind=4), dimension(dim), intent(in) :: phi
        real(kind=4), dimension(nx*nxx,ng), intent(in) :: SFPC,SF
        real(kind=4), dimension(nx*nxx) :: PF
        real(kind=4) :: som,val
        integer(kind=4) :: i,j,n
        open (10,file='app/Output/FLUX_SN.H')
        open (11,file='app/Output/PF_SN.H')
        call PowerPF(totNFM,Nmat,ng,nx,nxx,napc,SFPC,SF,PF)
        n=1;val=0.
        write(11,*) val
        do i = 1,nx
             do j=1,nxx
                write(11,*) PF(n)
                n=n+1  
             enddo
        enddo
        som = Delta(1)
        do i=2,totNFM
        write(10,*) som,(phi(i+j-1),j=0,dim-1,totNFM)
        som = som + Delta(i)
        enddo
        close(10)
        close(11)
end subroutine Plot_flux
!BL18
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
       print*, 'Dimension| One dimensions (1D)                              '
       print*, 'Geometry | Slab                                             '  
       print*, '____________________________________________________________'
       print*, ''
    end subroutine title1
!BL19
    subroutine title2()
       write(*,FMT=*)'**********************************************************'
       write(*,FMT=*)'                         Finished                         '                             
       write(*,FMT=*)'**********************************************************'  
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
!>BL21
    subroutine ScalarFluxPinC(dim,totNFM,Nmat,ng,nx,nxx,npx,npc,na,nfmesh,delta,&
                              assembly,phi,sigF,fmmid,core,SFPC,SF)
       !> CALCULATION OF SCALAR FLUX IN EACH PIN CELL
       integer(kind=4), intent(in) :: dim,totNFM,Nmat,ng,nx,nxx,npx,npc,na
       integer(kind=4), dimension(npc,npx), intent(in) ::  nfmesh
       real(kind=4), dimension(totNFM), intent(in) :: delta
       integer(kind=4), dimension(na,nxx), intent(in) :: assembly
       real(kind=4), dimension(dim), intent(in) :: phi
       integer(kind=4), dimension(nx), intent(in) :: core
       real(kind=4), dimension(Nmat,ng), intent(in) :: sigF
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=4), dimension(nx*nxx,ng), intent(out) :: SFPC,SF
       real(kind=4), dimension(totNFM,ng) :: flux
       integer(kind=4) :: i,j,k,l,n,n1
       real(kind=4) :: som
       open (10,file='app/Output/SFPC.H')
       n=1
       do i=1,ng
          do j=1,totNFM
             flux(j,i)=phi(n)
             n=n+1
          enddo
       enddo
       SFPC=0;SF=0   
       do j=1,ng
          n=1;n1=1
       do i=1,nx
             do k =1,nxx  
                do l=1,sum(nfmesh(assembly(core(i),k),:))
                   SFPC(n,j) =  SFPC(n,j) + flux(n1,j)  
                   SF(n,j)   =  SF(n,j)   + sigF(fmmid(n1),j)*flux(n1,j)
                   n1=n1+1
                enddo
              n=n+1
             enddo
       enddo
       enddo
       do i=1,nx*nxx
          write(10,*) (SFPC(i,j),j=1,ng)
       enddo
       close(10)
    end subroutine ScalarFluxPinC
!>BL22
    subroutine PowerPF(totNFM,Nmat,ng,nx,nxx,napc,SFPC,SF,PF)
       !> CALCULATION OF POWER PEAKING FACTOR
       implicit none
       integer(kind=4), intent(in) :: totNFM,Nmat,ng,nx,nxx,napc
       real(kind=4), dimension(nx*nxx,ng), intent(in) :: SFPC,SF
       real(kind=4), dimension(nx*nxx), intent(out) :: PF
       real(kind=4), dimension(nx*nxx) :: PF0
       real(kind=4) :: pm
       integer(kind=4) :: i,j,n
       PF0 = 0.0
       n=1
       do i=1,nx
          do j =1,nxx
             PF0(n) =  sum(SF(n,:))
             n=n+1
          enddo
       enddo
       pm  = sum(PF0)/float(napc)
       PF  = PF0/pm
    end subroutine PowerPF
!BLOC15
    subroutine NormalizeFlux(dim,totNFM,Nmat,ng,sigF,fmmid,delta,phi) 
    implicit none
    ! Neutron scalar flux is normalized according to sum(V*NusigF*phi=1)
    integer(kind=4), intent(in) :: dim,totNFM,Nmat,ng
    real(kind=4), dimension(Nmat,ng), intent(in) :: sigF
    real(kind=4), dimension(totNFM), intent(in) :: delta
    integer(kind=4), dimension(totNFM), intent(in) :: fmmid
    real(kind=4), dimension(dim), intent(inout) :: phi
    real(kind=4), dimension(totNFM,ng) :: flux
    integer(kind=4) :: i,j,k,n
    real(kind=4) :: norme,a1,a2,a3
    ! Initialize local variables
    n=1
       do i=1,ng
          do j=1,totNFM
             flux(j,i)=phi(n)
             n=n+1
          enddo
       enddo
    ! Normalized source     
    do i = 1,totNFM
        a1 = sum(flux(i,:)*delta(i)*sigF(fmmid(i),:))
        norme = norme  + sqrt(a1*a1)
    enddo
    flux = sum(delta)*(flux/norme)
    n=1
       do i=1,ng
          do j=1,totNFM
             phi(n) = flux(j,i)
             n=n+1
          enddo
       enddo
    end subroutine NormalizeFlux


