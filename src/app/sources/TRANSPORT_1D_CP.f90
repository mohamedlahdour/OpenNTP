!BL1
    subroutine Matrix_matrix_I(matrix_I,dim)
       implicit none
       integer(kind=4), intent(in) :: dim
       real(kind=8), dimension(dim,dim), intent(out) :: matrix_I 
       integer(kind=4) :: i,j
       do j = 1,dim
          do i = 1,dim
             if (i == j) then
             matrix_I(i,j) = 1.0D0
             else
             matrix_I(i,j) = 0.0D0
             endif
          enddo   
       enddo 
    end subroutine Matrix_matrix_I
!BL2
    subroutine fmm_id(nfmesh,RegMat,fmmid,nregion,totNFM)
       implicit none
       integer(kind=4), intent(in) :: nregion,totNFM
       integer(kind=4), dimension(nregion), intent(in) :: RegMat,nfmesh
       integer(kind=8), dimension(totNFM), intent(out) :: fmmid 
       integer(kind=4) :: m=0,i 
       do i = 1,nregion
             fmmid(m+1:m+nfmesh(i)) = RegMat(i)
             m = m+ nfmesh(i)
       enddo  
    end subroutine fmm_id
!BL3
    function    en(n,xx) 
!      Exponential integrals E1 through E4
!      Range    -- all positive arguments.
!      Accuracy -- absolute accuracy E-06 to E-07.
!      Method   -- for E1 two chebyshev expansions are used, one below 
!                  and the other above 4.
!                  for E2 through E4 the following recurrence relation is used:
!                  EN(X) = (EN(-X)-X*E(N-1)(X))/(N-1)
       implicit none
       real(kind=8), intent(in) :: xx
       integer(kind=4), intent(in) :: n
       real(kind=8) :: t,f,r,s,ex,enn=0.
       real(kind=8) :: en
       real(kind=8) :: x
       integer(kind=4) :: i
       real(kind=8), dimension(6) :: A
       real(kind=8), dimension(4) :: B
       real(kind=8), dimension(4) :: C
       data  A / -0.5772156649,0.99991930,-0.24991055,0.05519968,&
                 -0.0097600400,0.00107857/
       data  B /  8.5733287401,18.0590169730,8.63476082500,0.2677737343/
       data  C /  9.5733223400,25.6329561486,21.0996530827,3.9584969228/
       x = abs(xx) + 1.0E-20
       if ( n < 1 ) then
       write (*,'(a)') 'error (n in function EN(n,x) must be >= 0)'
       elseif ( x < -1.0E-10 ) then
       write (*,'(a)') 'error (x in function EN(n,x) must be >= 0)'
       elseif ( x <= 0 ) then
           if ( n <= 1 ) then
           enn = 1.E20
           else
           enn = 1/(n-1)
           endif
       else
       ex=exp(-x)
           if ( n == 0 ) then
           enn = ex/x
           else
              if ( x <= 1 ) then
              t  = a(1)+x*(a(2)+x*(a(3)+x*(a(4)+x*(a(5)+x*a(6)))))
              f  = t-log(x)
              enn = f
              else
              r  = b(4)+x*(b(3)+x*(b(2)+x*(b(1)+x)))
              s  = c(4)+x*(c(3)+x*(c(2)+x*(c(1)+x)))
              f  = r/s*ex/x
              enn  = f
              endif
              do  i=1,n-1
              f =(ex-x*f)/i
              end do
              enn = f
           endif
       endif
       en = enn
    end function en
!BL4   
    function fii(signe,taux,SigTi,li)
       implicit none
       real(kind=8), intent(in) :: signe
       real(kind=8), intent(in) :: SigTi, li, taux
       real(kind=8) :: fii,a0,a1,a2,a3,EN
       if     ( SigTi /= 0) then
       a0   = 1/(li*SigTi*SigTi) 
       a1   = EN(2,taux)   
       a2   = EN(3,taux) 
       a3   = EN(3,taux +signe*SigTi*li) 
       fii  = (signe*a1)/SigTi-a0*(a2-a3)
       else
       fii  = 0.5*li*EN(1,taux)
       end if
    end function fii
!BL5
    function fij(signe,taux,li,lj,SigTi,SigTj)
       implicit none
       real(kind=8), intent(in) :: signe
       real(kind=8), intent(in) :: li,lj,SigTi,SigTj,taux
       real(kind=8) :: fij,a0,a1,a2,a3,a4,EN
       if     ( SigTi /= 0 .and. SigTj /= 0 ) then
       a0   = 1/(li*SigTi*SigTj)
       a1   = EN(3,taux)
       a2   = EN(3,taux+signe*li*SigTi)
       a3   = EN(3,taux+signe*lj*SigTj)
       a4   = EN(3,taux+signe*li*SigTi+signe*lj*SigTj)
       fij  = 0.5*a0*(a1-a2-a3+a4)
       elseif ( SigTi == 0 .and. SigTj /= 0 ) then
       a0   = 1/SigTj
       a1   = EN(2,taux)
       a2   = EN(2,taux+signe*lj*SigTj)
       fij  = signe*0.5*li*a0*(a1-a2)
       elseif ( SigTi /= 0 .and. SigTj == 0 ) then
       a0   = 1/SigTi
       a1   = EN(2,taux)
       a2   = EN(2,taux+signe*li*SigTi)
       fij  = signe*0.5*lj*a0*(a1-a2)
       else
       fij  = 0.5*li*lj*EN(1,taux)
       end if
    end function fij
!BL6
    subroutine Transport_Corr(SigSS,SigTT,SigS,SigT,ng,Nmat,order)
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,order
       real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigSS
       real(kind=8), dimension(Nmat,ng), intent(in) :: SigTT
       real(kind=8), dimension(Nmat,ng,ng), intent(out) :: SigS
       real(kind=8), dimension(Nmat,ng), intent(out) :: SigT
       integer(kind=4) :: i,j
       if (order==1) then
          SigT = SigTT
          SigS(:,:,:) = SigSS(:,1,:,:)
       else
          do i = 1,Nmat
             do  j = 1,ng
             SigT(i,j) = SigTT(i,j) - SigSS(i,2,j,j)
             SigS(i,j,j) = SigSS(i,1,j,j) - SigSS(i,2,j,j)
             enddo
          enddo
       endif
    end subroutine Transport_Corr
!BL7
    subroutine  Pij_f(vol,albedo,SigT,fmmid,pij,ng,dim,totNFM,Nmat)
       integer(kind=4), intent(in) :: ng,dim,totNFM,Nmat
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid 
       real(kind=8), dimension(dim,dim), intent(out) :: Pij
       real(kind=8), dimension(totNFM), intent(in) :: vol
       real(kind=8), dimension(Nmat,ng), intent(in) :: SigT
       real(kind=8), dimension(2), intent(in) :: albedo

       real(kind=8), dimension(dim) :: sig,vo
       real(kind=8), dimension(dim,dim) :: test,pij1
       real(kind=8) :: a0,a1,taucell,tau0,fii,fij
       integer(kind=4) :: ri,rj,k,k0=1,k1,K2,i,j
       pij  = 0.0D0
       pij1 = 0.0D0
       i = 1
       do k1=1,ng
       do K2=1,totNFM
       sig(i) = SigT(fmmid(k2),k1)
       vo(i)  = vol(k2)
       i = i + 1
       enddo
       enddo

       if  ( albedo(1) == 0.0D0 .and. albedo(2) == 0.0D0 ) then
           do k = 1, ng
              do  ri = k0,totNFM*k
                  pij(ri,ri) = fii(1.0D0,0.0D0,sig(ri),vo(ri))
                  tau0 = 0.0
                  do  rj = ri+1,totNFM*k
                      pij(ri,rj) = fij(1.0D0,tau0,vo(ri),vo(rj),sig(ri),sig(rj))
                      tau0 = tau0 + vo(rj)*sig(rj)
                  enddo
              enddo
              k0 = totNFM+k0
           enddo

       elseif ( albedo(1) == 1.0D0 .and. albedo(2) == 1.0D0 ) then
           k0 = 1
       do  k = 1,ng
           m = - 1
           taucell = DOT_PRODUCT(SigT(fmmid(:),k),vol) 
           test = 1.1
           do while ( maxval(test(:,:)) >= 1.0E-10 ) 
                    m   = m + 1 
                        tau0 = 0.0
                    do  ri = k0,totNFM*k
                        a0 = albedo(1)**(m)*fii(1.0D0,m*taucell,sig(ri),vo(ri))
                        a1 = albedo(2)**(m+1)*fii(-1.0D0,(m+1)*taucell,sig(ri),vo(ri))
                        pij1(ri,ri) = pij1(ri,ri) + a0 + a1
                        tau0 = 0.0
                        do  rj = ri+1,totNFM*k
                        a0 = albedo(1)**(m)*fij(1.0D0,(m*taucell+tau0),vo(ri),vo(rj),sig(ri),sig(rj))
                        a1 = albedo(2)**(m+1)*fij(-1.0D0,((m+1)*taucell-tau0),vo(ri),vo(rj),sig(ri),sig(rj))
                            pij1(ri,rj) = pij1(ri,rj) + a0 + a1
                            tau0 = tau0 + vo(rj)*sig(rj) 
                        enddo
                    enddo
           pij  = pij + pij1 
           test = pij1
           pij1 = 0.
           enddo
       k0 = totNFM+k0
       enddo
       else
       print*, "This limit condition is not available"
       stop
       endif
       k0=1
       do k = 1,ng
       do i = k0,totNFM*k
          do j = i,totNFM*k  
             Pij(j,i) = (vo(i)/vo(j))*Pij(i,j)
          end do
       end do 
       k0 = totNFM  + k0
       enddo
    end subroutine Pij_f
!BL8
    subroutine flux_guess(NusigF,vol,fmmid,phi_guess,dim,ng,totNFM,Nmat)
       implicit none
       integer(kind=4), intent(in) :: ng,totNFM,Nmat,dim
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid 
       real(kind=8), dimension(Nmat,ng), intent(in) :: NusigF
       real(kind=8), dimension(totNFM), intent(in) :: vol
       real(kind=8), dimension(dim), intent(out) :: phi_guess
       real(kind=8), dimension(dim) :: a10,a11
       real(kind=8) :: kguess
       integer(kind=4) :: i,m,n
       kguess = 1.0
       i = 1
       do m=1,ng
          do n=1,totNFM
          a10(i) = NusigF(fmmid(n),m)
          a11(i) = vol(n)
          i = i + 1
          enddo
       enddo 
       phi_guess  = kguess/dot_product(a10,a11)
    end subroutine flux_guess
!BL9
    subroutine Matrix_D(SigS,fmmid,D,ng,totNFM,Nmat,dim)
       implicit none
       integer(kind=4), intent(in) :: ng,totNFM,Nmat,dim
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid 
       real(kind=8), dimension(Nmat,ng,ng), intent(in) :: SigS
       real(kind=8), dimension(dim,dim), intent(out) :: D
       integer(kind=4) :: k0,k,i,j
            k0 = 1
            do   k   = 1, ng
                 j   = 1
                 do  i  = k0,totNFM*k
                     D(i,i) = SigS(fmmid(j),k,k)
                     j = j + 1
                 end do
            k0 = totNFM + k0
            enddo
    end subroutine Matrix_D
!BL10
    subroutine Matrix_C(matrix_I,D,pij,C,dim)
       implicit none
       integer(kind=4), intent(in) :: dim
       real(kind=8), dimension(dim,dim), intent(in) :: matrix_I,D,pij
       real(kind=8), dimension(dim,dim), intent(out) :: C
       C(:,:) = matrix_I(:,:) - matmul(D,pij)
    end subroutine Matrix_C 
!BL11
    subroutine matinv(a,ainv,n)
       implicit none
       integer(kind=4), intent(in) :: n
       real(kind=8), dimension(n,n), intent(in) :: a
       real(kind=8), dimension(n,n), intent(out) :: ainv
       real(kind=8), dimension(n,2*n) :: b
       integer(kind=4) :: i,j,k
       real(kind=8) :: pivot=0.0D0,xnum
!      make augmented matrix 
       do i=1,n 
       do j=1,n 
       b(i,j)=0.0
       b(i,j+n)=0.0
       b(i,j)=a(i,j)
       if(i.eq.j) then 
       b(i,j+n)=1.0
       end if  
       end do
       end do 
       do i=1,n
!      choose the leftmost non-zero element as pivot 
       do j=1,n
       if (abs(b(i,j)).gt. 0.0) then
       pivot=b(i,j)
       exit 
       end if
       end do 
!      step 1: change the chosen pivot into "1" by dividing 
!      the pivot's row by the pivot number 
       do j=1,2*n
       b(i,j)=b(i,j)/pivot
       end do
       pivot=b(i,i)  !update pivot value 
!      step 2: change the remainder of the pivot's column into 0's
!      by adding to each row a suitable     multiple of the pivot row 
       do k=1,n !row 
       if(k.ne.i) then
       xnum=b(k,i)/pivot   !same column with the current pivot
       do j=1,2*n !col 
       b(k,j)=b(k,j)-xnum*b(i,j) 
       end do 
       end if 
       end do 
       end do 
!      prepare the final inverted matrix 
       do i=1,n 
       do j=1,n 
       ainv(i,j)=b(i,j+n) 
       end do 
       end do 
       return
    end subroutine matinv 
!BL12
    subroutine Matrix_W(ainv,pij,W,dim)
       implicit none
       integer(kind=4), intent(in) :: dim
       real(kind=8), dimension(dim,dim), intent(in) :: ainv,pij
       real(kind=8), dimension(dim,dim), intent(out) :: W
       W(:,:) = matmul(pij,ainv)
    end subroutine Matrix_W
!BL13
    subroutine Matrix_L(SigS,fmmid,L,ng,totNFM,Nmat,dim)
       implicit none
       integer(kind=4), intent(in) :: ng,totNFM,Nmat,dim
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid 
       real(kind=8), dimension(Nmat,ng,ng), intent(in) :: SigS
       real(kind=8), dimension(dim,dim), intent(out) :: L
       integer(kind=4) :: i,k0,k1,k2,k3
                   L(:,:) = 0.0
                k2  = 0
            do  k1  = 1,ng
                k0  = 1
                k3  = 0
                    do while (k0<k1) 
                             do  i  = 1,totNFM
                             L(i+k3,i+(k1-1)*totNFM) = SigS(fmmid(i),k0,k1) ! 1 --> 2 SigS(1,1,2)
                             enddo
                             k0 = k0 + 1
                             k3 = k3 + totNFM
                    enddo
            enddo
    end subroutine 
!BL14    
    subroutine Matrix_U(SigS,fmmid,U,ng,totNFM,Nmat,dim)
       implicit none
       integer(kind=4), intent(in) :: ng,totNFM,Nmat,dim
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid 
       real(kind=8), dimension(Nmat,ng,ng), intent(in) :: SigS
       real(kind=8), dimension(dim,dim), intent(out) :: U
       integer(kind=4) :: i,k0,k1,k2,k3
                U(:,:) = 0.0
                k2  = 0
            do  k1  = 1,ng
                k0  = 1
                k3  = 0
                    do while (k0<k1) 
                             do  i  = 1,totNFM
                             U(i+(k1-1)*totNFM,i+k3) = SigS(fmmid(i),k1,k0) ! 1 --> 2 SigS(1,1,2)
                             enddo
                             k0 = k0 + 1
                             k3 = k3 + totNFM
                    enddo
            enddo
    end subroutine Matrix_U
!BL15
    subroutine Matrix_F(NusigF,CHI,fmmid,F,ng,totNFM,Nmat,dim)
       implicit none
       integer(kind=4), intent(in) :: ng,totNFM,Nmat,dim
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid 
       real(kind=8), dimension(Nmat,ng), intent(in) :: NusigF,CHI
       real(kind=8), dimension(dim,dim), intent(out) :: F
       integer(kind=4) :: i,j,k0,k1,k2,k3,k
            F(:,:) = 0.0
            k0 = 1
            do   k   = 1, ng
                 j   = 1
                 do  i  = k0,totNFM*k
                     F(i,i) = CHI(fmmid(j),k)*NusigF(fmmid(j),k)
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
                             F(i+k3,i+(k1-1)*totNFM) =  CHI(fmmid(i),k1)*NusigF(fmmid(i),k0)
                             F(i+(k1-1)*totNFM,i+k3) =  CHI(fmmid(i),k0)*NusigF(fmmid(i),k1)
                             enddo
                             k0 = k0 + 1
                             k3 = k3 + totNFM
                    enddo
            enddo
    end subroutine Matrix_F
!BL16
    subroutine Matrix_AB(matrix_I,L,W,U,F,A,B,dim)
       implicit none
       integer(kind=4), intent(in) :: dim
       real(kind=8), dimension(dim,dim), intent(in) :: matrix_I,L,W,U,F
       real(kind=8), dimension(dim,dim), intent(out) :: A,B
       A(:,:)   = matrix_I(:,:) - matmul(W,L) - matmul(W,U)
       B(:,:)   = matmul(W,F)
    end subroutine Matrix_AB
!BL17
    subroutine aleig(a,b,eps,iter,eval,phi,phi_guess,ng,totNFM,Max_it,dim)
!      -------------------  The inverse power method  -------------------
!                           (A - 1/Keff*B)*phi = 0
! 
!      ------------------------------------------------------------------
       implicit none
       integer(kind=4), intent(in) :: ng,totNFM,Max_it,dim
       real(kind=8), intent(in) :: eps
       real(kind=8), dimension(dim,dim), intent(in) :: a,b
       real(kind=8), dimension(dim), intent(in) :: phi_guess
       real(kind=8), dimension(dim), intent(out) :: phi
       real(kind=8), intent(out) :: eval
       integer(kind=4), intent(out) :: iter
       real(kind=8), dimension(dim,dim) :: ai,ainv
       real(kind=8), dimension(dim) :: gar,vec1,vec2
       real(kind=8), dimension(dim) :: phi1,phi2
       real(kind=8), dimension(ng) :: moy
       real(kind=8) :: s1,s2,err1,err2,test,eval1,eval2,epsil1=1.0,epsil2=1.0,norme
       integer(kind=4) :: i,k0,k1
       call matinv(a,ainv,dim)
       ai = matmul(ainv,b)
       iter = 0
       phi1 = phi_guess
       eval1 = 0.0
       do while ( epsil1 >= eps .and. epsil2 >= eps )       
           iter = iter + 1
           if ( iter > Max_it ) then
           print*, 'error(unable to converge(1))'
           stop
           end if
           gar    = matmul(ai,phi1)
           vec1   = matmul(a,phi1)
           vec2   = matmul(b,phi1)
           s1     = dot_product(vec1,vec2)
           s2     = dot_product(vec2,vec2)
           eval2  = s1/s2
!          critere de convergence sur les valeurs porpres               
           epsil1 =  abs(eval2-eval1)/eval2
           eval1  = eval2 
           phi2   = gar*eval1         
           err1   = maxval(abs(phi2))
           err2   = maxval(abs(phi2-phi1))
!          critere de convergence sur les vecteurs porpres   
           epsil2 = err2/err1       
           phi1   = phi2      
           if (iter == 1)  then
                   test = epsil1
           elseif (iter >= 10 .and. epsil1 > test) then
           print*, 'error(unable to converge(2) )'
           print*, 'because error in iteration ten is sup to one'
           end if
           norme = sqrt(sum(phi1*phi1))
           phi1 = phi1/norme
           phi = phi1
           eval = eval1
           write(*,2000)iter,1/eval,epsil1
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
    end subroutine aleig
!BL18
    subroutine vol_ray(entree,nfmesh,rayon,ray,vol,nregion,totNFM)
        implicit none
        integer(kind=4), intent(in) :: totNFM,nregion
        integer(kind=4), dimension(nregion), intent(in) :: nfmesh
        real(kind=8), dimension(nregion), intent(in) :: rayon
        real(kind=8), dimension(totNFM), intent(out) :: ray,vol
        integer(kind=4) :: k1,k2=0,i
        real(kind=8) :: som=0.,PI,cte = 4.0d0/3.0d0
        logical entree
        data PI     /3.1415927/    
        do  k1 = 1,nregion
            if (entree .eqv. .false.) then
               vol(k2+1:nfmesh(k1)+k2) = (PI*rayon(k1)**2-som)/float(nfmesh(k1))
               som = PI*rayon(k1)**2
            else
               vol(k2+1:nfmesh(k1)+k2) = (cte*PI*rayon(k1)**3-som)/float(nfmesh(k1))
               som = cte*PI*rayon(k1)**3
            endif 
 
            if (entree .eqv. .false.) then
                do   i = 1,nfmesh(k1)
                     ray(k2+i) = sqrt(sum(vol(1:i+k2))/PI)
                enddo
            else
                do   i = 1,nfmesh(k1)
                     ray(k2+i) = (sum(vol(1:i+k2))/(cte*PI))**(1.0d0/3.0d0)
                enddo
            endif
            k2 = k2 + nfmesh(k1)
        enddo

    end subroutine vol_ray
!BL19
       function kin(n,gho)
       implicit none
       real(kind=8), intent(in) :: gho
       integer(kind=4), intent(in) :: n
!      first and second-order bickley functions  ki1 et ki2
!      ki1(x) interval boundaries: 0.04 - 6.0 - 15.0; A() and C() coeffs
!      ki2(x) interval boundaries: 0.50 - 4.0 - 15.0; B()  coeffs
!      ki3(x)  3 and 4-term power series corresponding to chebyshev expansions
!      accuracy 1.0E-5 up to argument 9.
!      intervals of 0.05 in range from 0.0 to 0.5, coeffs D();
!                   0.10 in range from 0.5 to 1.0, coeffs D();
!                   0.40 in range from 1.0 to 5.0, coeffs E();
!                   0.80 in range from 5.0 to 9.0, coeffs E().
       real(kind=8), dimension(27) :: A
       real(kind=8), dimension(45) :: B
       real(kind=8),  dimension(3) :: C
       real(kind=8), dimension(45) :: D
       real(kind=8), dimension(51) :: E
       integer(kind=8), dimension(20) :: INDEXD,INDEXE
       integer(kind=8) :: l,i,imin
       real(kind=8) :: u=0.0D0,f1,a0=0.0D0,a1,a2=0.0D0,kin,x2,pi=4*atan(1.0)
       real(kind=8) :: b0,b1,b2
       real(kind=8) :: x
       data A /10.584648, 6.5454700,1.97645700, 0.3334960,&
                0.035407, 0.0025680,0.00013500, 0.0000050,&
               10.762666, 5.6233350,1.43543700, 0.2125040,&                               
                0.020365, 0.0013600,0.00006700, 0.0000030,&
               0.0004364,-0.0003873,0.00027400,-0.0001581,&
               0.0000763,-0.0000315,0.00001130,-0.0000036,&
               0.0000010,-0.0000003,0.00000010/
       data B /1.4534664,-0.2436620,0.02584650,-0.0029653,&
               0.0005322,-0.0001499,0.00005600,-0.0000249,&
               0.0000125,-0.0000068,0.00000400,-0.0000025,&
               0.0000017,-0.0000012,0.00000090,-0.0000008,&
               0.0000007, 0.3039967,-0.2136079, 0.0961280,&
              -0.0324165, 0.0091054,-0.0023228, 0.0005813,&
              -0.0001516, 0.0000426,-0.0000129, 0.0000042,&
              -0.0000014, 0.0000005,-0.0000002, 0.0000001,&
               0.0031373,-0.0028564,0.00216700,-0.0013874,&
               0.0007615,-0.0003642,0.00015410,-0.0000585,&
               0.0000202,-0.0000064,0.00000190,-0.0000005,&
               0.0000001/
       data C /0.4227843, 0.0534107,0.08333333/
       data D /0.7853961,-0.9990226, 0.7266088, 0.7852024,&
              -0.9912340, 0.6466375, 0.7845986,-0.9791293,&
               0.5856605, 0.7834577,-0.9638914, 0.5346648,&
               0.7817094,-0.9463843, 0.4907827, 0.7793031,&
              -0.9271152, 0.4521752, 0.7762107,-0.9054822,&
               0.4177388, 0.7724519,-0.8849865, 0.3869945,&
               0.7679903,-0.8626685, 0.3590753, 0.7628988,&
              -0.8400133, 0.3338676, 0.7540982,-0.8054172,&
               0.2998569, 0.7401279,-0.7587821, 0.2609154,&
               0.7239594,-0.7125290, 0.2278226, 0.7058777,&
              -0.6672761, 0.1994999, 0.6861762,-0.6234536,&
               0.1751248/
       data E /0.724729400,-0.753835500, 0.320322300,-5.337485e-2,&
               0.666372000,-0.627975200, 0.229528000,-3.146833e-2,&
               0.595616300,-0.509412400, 0.163166700,-1.906198e-2,&
               0.519103100,-0.404600700, 0.115241800,-1.174720e-2,&
               0.442595400,-0.315954800, 8.097913e-2,-7.328415e-3,&
               0.370317800,-0.243434100, 5.669960e-2,-4.617254e-3,&
               0.168402200,-7.158569e-2, 7.923547e-3, 0.127830700,&
              -5.016344e-2, 5.095111e-3, 9.611422e-2,-3.501524e-2,&
               3.286040e-3, 7.170491e-2,-2.437465e-2, 2.126242e-3,&
               4.616317e-2,-1.425519e-2, 1.123687e-3, 2.475115e-2,&
              -6.810124e-3, 4.762937e-4, 1.302864e-2,-3.232035e-3,&
               2.031843e-4, 6.742673e-3,-1.524126e-3, 8.701440e-5,&
               3.454768e-3,-7.157367e-4, 3.742673e-5/
       data INDEXD / 3,  6,  9, 12, 15, 18, 21, 24, 27, 30,&
                    33, 33, 36, 36, 39, 39, 42, 42, 45, 45/
       data INDEXE / 4,  8, 12, 16, 20, 24, 27, 30, 33, 36,&
                    39, 39, 42, 42, 45, 45, 48, 48, 51, 51/

!      ----------------------- test sur gho et n ----------------------  
       if (n < 1) then 
       write (*,'(a)') 'error (n must be > 0)'
       stop
       elseif (gho < -1.e-10) then
       write (*,'(a)') 'error (gho must be >= 0)'
       stop
       endif
!      ----------------------------------------------------------------
       x = abs(gho)
       go to (10,90,91)N
!      ---------------------------- ki1(x) ----------------------------
       10  if (x .gt. 1.e-5) go to 20
           kin  = pi/2
           return
       20  if (x .ge. 0.04) go to 30
           x2   = x*x
           kin  = pi/2-x*( (x2*C(2)*C(1))-alog(0.5*real(x,4))*(x2*C(3)+1.0) )
           return
       30  if (x .ge. 6.00) go to 50
           u    = 0.11111111*x*x-2.0
           l    = 1
           i    = 8
           imin = 1
           go to 70
       40  f1   = 0.5*(a0-a2)
           l    = 2
           i    = 16
           imin = 9
           go to 70
       50  if (x .gt. 15.0) go to 60
           u    = 0.44444444*x-4.6666667
           l    = 3
           i    = 27
           imin = 17
           go to 70
       60  kin  = 0.0
           return
       70  a0   = A(i)
           a1   = 0.0
           a2   = 0.0
       80  i    = i-1
           a2   = a1
           a1   = a0
           a0   = u*a1-a2+a(i)
           if (i .gt. imin) go to 80
           if (l .eq. 1) go to 40
           kin  = 0.5*(a0-a2)
           if (l .eq. 2) kin=pi/2-x*(f1-kin*alog(0.5*real(x,4)))
           return
!      ---------------------------- ki2(x) ----------------------------
       90  if (x .ge. 0.5) go to 100
           u    = 8.0*x-2.0
           i    = 17
           imin = 1
           go to 120
       100 if (x .ge. 4.0) go to 110
           u    = 1.1428571*x-2.5714286
           i    = 32
           imin = 18
           go to 120
       110 if (x .gt. 15.0) go to 40
           u    = 0.36363636*x-3.4545455
           i    = 45
           imin = 33
       120 b0   = B(i)
           b1   = 0.0
           b2   = 0.0
       130 i    = i-1
           b2   = b1
           b1   = b0
           b0   = u*b1-b2+B(i)
           if (i .gt. imin) go to 130
           kin  = 0.5*(b0-b2)
           return
!      ---------------------------- ki3(x) ----------------------------
       91  if (x .ge. 1.0) go to 140
           i    = ifix(20.0*real(x,4))+1
           i    = INDEXD(i)
           kin  = x*(x*D(i)+D(i-1))+D(i-2) 
           return
       140 i    = ifix(2.5*(real(x,4)-1.0))+1
           if (x .ge. 3.4) go to 150
           i    = INDEXE(i)
           kin  = x*(x*(x*E(I)+E(i-1))+E(i-2))+E(i-3)
           return
       150 if (x .ge. 9.0) go to 160
           i    = INDEXE(i)
           kin  = x*(x*E(i)+E(i-1))+E(i-2)
           return
       160 kin  = 0.0
           return
    end function kin
!BL20
    function Di_f(sigti,li,entree)
       implicit none
       real(kind=8) :: Di_f,kin,pi=4*atan(1.0)
       real(kind=8), intent(in) :: li
       real(kind=8), intent(in) :: sigti
       logical entree
       if (entree .eqv. .false.) then
           if    (sigti .ne. 0.0)  then
           Di_f  = li/sigti-(kin(3,0.0D0)-kin(3,sigti*li))/sigti**2
           else
           Di_f  = (li**2*pi)/4
           end if
       else
           if    (sigti .ne. 0.0)  then
           Di_f  = li/sigti-(1-exp(-sigti*li))/sigti**2
           else
           Di_f  = 0.5*li**2
           end if
       end if
    end function Di_f     
!BL21
    function Cij_f(gho,sigti,sigtj,li,lj,entree)   
       implicit none
       real(kind=8) :: Cij_f,kin
       real(kind=8), intent(in) :: gho
       real(kind=8), intent(in) :: li,lj
       real(kind=8), intent(in) :: sigti,sigtj
       logical entree
       if (entree .eqv. .false.) then
            if     (sigti .ne. 0.0 .and. sigtj .ne. 0.0) then
            Cij_f = (kin(3,gho)-kin(3,gho +sigti*li)-&
               kin(3,gho +sigtj*lj)+&
               kin(3,gho +sigti*li+sigtj*lj))/(sigti*sigtj)
            elseif (sigti .eq. 0.0 .and. sigtj .ne. 0.0) then
            Cij_f = (kin(2,gho)-&
               kin(2,gho +sigtj*lj))*(li/sigtj)
            elseif (sigti .ne. 0.0 .and. sigtj .eq. 0.0) then
            Cij_f = (kin(2,gho)-&
               kin(2,gho +sigti*li))*(lj/sigti)
            else   
            Cij_f = li*lj*kin(1,gho)
            end if
        else
            if     (sigti .ne. 0.0 .and. sigtj .ne. 0.0) then
            Cij_f = (exp(-gho)-exp(-(gho +sigti*li))-&
               exp(-(gho +sigtj*lj))+&
               exp(-(gho +sigti*li+sigtj*lj)))/(sigti*sigtj)
            elseif (sigti .eq. 0.0 .and. sigtj .ne. 0.0) then
            Cij_f = (exp(-gho)-&
               exp(-(gho +sigtj*lj)))*(li/sigtj)
            elseif (sigti .ne. 0.0 .and. sigtj .eq. 0.0) then
            Cij_f = (exp(-gho)-&
               exp(-(gho +sigti*li)))*(lj/sigti)
            else   
            Cij_f = li*lj*exp(-gho)
            end if
        endif
    end function Cij_f
!BL22
    subroutine tracking_f(entree,ray,track1,track3,totNFM)
       implicit none
       integer(kind=4), intent(in) :: totNFM
       real(kind=8), dimension(totNFM), intent(in) :: ray
       real(kind=8), dimension(totNFM,totNFM,2), intent(out) :: track3
       real(kind=8), dimension(totNFM,2), intent(out) :: track1
       real(kind=8), dimension(2) :: x,w
       real(kind=8), dimension(totNFM) :: z
       real(kind=8) :: rik1,rik2,r,rd,aux,ct1,ct2,pi
       integer(kind=4) :: i,ik,il,m,ngauss=2
       logical entree
       data pi /3.1415927/
       data x /0.3550510257,0.8449489743/
       data w /0.1819586183,0.3180413817/
       !data x /0.0985350858, 0.3045357266, 0.5620251898, 0.8019865821, 0.9601901429/
       !data w /0.0157479145, 0.0739088701, 0.1463869871, 0.1671746381, 0.0967815902/
       !data x /.0730543287,.2307661380,.4413284812,.6630153097,.8519214003,.9706835728/
       !data w /.0087383018,.0439551656,.0986611509,.1407925538,.1355424972,.0723103307/
       rik1=0.0
       do ik = 1,totNFM    
          rik2 = ray(ik)  
          rd = rik2-rik1  
          do il = 1,ngauss
             r = rik2-rd*x(il)**2
             aux = 4.*rd*w(il)
             if (entree .eqv. .true.) then
             aux = aux*r*pi
             endif
             ct1=0.0
             m=totNFM -ik+1
             do i = ik,totNFM 
                ct2 = sqrt(ray(i)**2-r**2)
                z(i-ik+1) =ct2-ct1
                ct1 = ct2
             enddo
             track3(ik,:,il) = z(:)
             track1(ik,il) = aux                              
          end do
          rik1=rik2
       end do
       end subroutine tracking_f
!BL23
   subroutine  Pij_Cyl_Sph(entree,SigT,fmmid,Pij_table,track1,track3,&
                                           vol,ray,albedo,ng,totNFM,Nmat,dim)
       implicit none
       integer(kind=4), intent(in) :: albedo,ng,totNFM,Nmat,dim
       real(kind=8), dimension(dim,dim), intent(out) :: Pij_table
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=8), dimension(Nmat,ng), intent(in) :: SigT
       real(kind=8), dimension(totNFM,totNFM,2),intent(in) :: track3
       real(kind=8), dimension(totNFM,2),intent(in) :: track1
       real(kind=8), dimension(totNFM),intent(in) :: vol,ray
       real(kind=8), dimension(dim) :: sig,vo
       real(kind=8), dimension(totNFM) :: seg
       real(kind=8), dimension(dim,dim) :: Pij
       real(kind=8), dimension(dim) :: Pis,psi,ones
       real(kind=8), dimension(ng) :: Pss,beta
       real(kind=8) :: pi,val,pwr,tau,tau1j,tauij,Di_f,Cij_f,sum
       integer(kind=4):: ix,i,j,ip,jp,k0=1,k,k1,k2,ir
       logical entree
       Pij_table = 0.d0
       Pij = 0.d0
       data pi     /3.1415927/
       data val    /0.0/
       if (albedo == 2) then
       print*, "This limit condition is not available"
       stop
       else
       ones = 1.0
       i = 1
       do k1=1,ng
          do K2=1,totNFM
             sig(i) =  SigT(fmmid(k2),k1)
             vo(i)  = vol(k2)
             i = i + 1
          enddo
       enddo 

       do k=1,ng
          ir = 1
          do  ix = k0,totNFM*k
              do  i = 1,2                ! On somme sur le 1ere track
                  seg  = track3(ir,:,i)   
                  pwr  = track1(ir,i)
                  Pij(ix,ix) = Pij(ix,ix) + pwr*Di_f(sig(ix),2.0D0*seg(1),entree)       
                  tau=2.0D0*seg(1)*sig(ix)
                  tau1j=0.0            
                  do  ip = ix+1,totNFM*k           ! On somme sur le 2eme track 
                      Pij(ip,ip) = Pij(ip,ip) +pwr*(2.0D0*Di_f(sig(ip),& 
                                   seg(ip-ix+1),entree)+ &
                                   Cij_f(tau,sig(ip),sig(ip), &
                                   seg(ip-ix+1),seg(ip-ix+1),entree)) 
                      Pij(ix,ip) = Pij(ix,ip) + &
                                   pwr*(Cij_f(tau1j,sig(ix),sig(ip),&
                                   2.0D0*seg(1),seg(ip-ix+1),entree)) 
                                   tauij=0.0D0                      
                      do  jp = ip+1,totNFM*k
                          Pij(ip,jp) = Pij(ip,jp) + pwr*(Cij_f((tauij+tau+seg(ip-ix+1)* &
                                       sig(ip)),sig(ip),sig(jp),seg(ip-ix+1), &
                                       seg(jp-ix+1),entree)+ &
                                       Cij_f(tauij,sig(ip),sig(jp),&
                                       seg(ip-ix+1),seg(jp-ix+1),entree))
                          tauij = tauij + seg(jp-ix+1)*sig(jp)
                      end do
                      tau  = tau + 2.0D0*seg(ip-ix+1)*sig(ip) ;
                      tau1j = tau1j+ seg(ip-ix+1)*sig(ip)
                  
                  end do
             end do
          ir = ir + 1 
          end do         
    
       do i = k0,totNFM*k
          do j = i,totNFM*k
             val = Pij(i,j)
             Pij(i,j) = val/vo(i)
             Pij(j,i) = val/vo(j)
          end do
       end do 
       k0 = totNFM+k0
       end do
       
!      Pis non reduite 
       Pis = ones - matmul(Pij,sig)        
!      psi reduite                    
       psi =  (4.0D0*vo*Pis/(2.0D0*pi*ray(totNFM)))                
!      Pss non reduite       
       k0 = 1
       do k2 = 1,ng
          sum = 0.0D0  
          do j = k0,totNFM*k2
             sum = sum + psi(j)*sig(j)
          enddo
          Pss(k2) =  1.0D0 - sum
          k0 = totNFM+k0
       enddo    
!      Pij_table assume a white boundary condition
       do k2 = 1,ng
       beta(k2) = albedo/(1.0D0-albedo*Pss(k2)) 
       enddo  

       k0 = 1
       do k = 1,ng
       do i = k0,totNFM*k
          do j = k0,totNFM*k
               Pij_table(i,j) = Pij(i,j) + beta(k)*Pis(i)*psi(j)
          end do     
       end do
          k0 = totNFM+k0
       enddo   
       endif

    end subroutine  Pij_Cyl_Sph
!BL24
    subroutine plot_Sla_Cyl_Sph(vol,expression,phi,totNFM,dim)
        implicit none
        integer(kind=4), intent(in) :: dim,totNFM
        real(kind=8), dimension(dim), intent(in) :: phi
        real(kind=8), dimension(totNFM), intent(in) :: vol
        logical expression
        integer(kind=4) :: i,j
        real(kind=8) :: som
        som=vol(1)
        open (11,file='app/Output/flux_cp.h')
        if (expression .eqv. .true.) then
        do i=1,totNFM
        write(11,'(1p6e12.4)') som,(phi(i+j),j=0,dim-1,totNFM)
        som = som + vol(i)       
        enddo
        elseif (expression .eqv. .false.) then
        do i=1,totNFM
        write(11,'(1F10.6,200F10.6)') vol(i),(phi(i+j),j=0,dim-1,totNFM) 
        enddo
        endif
        close(11)
    end subroutine plot_Sla_Cyl_Sph
!BL25
subroutine Output(start,albedo,tm,k_eff,SigT,NusigF,SigS,Chi,dcell,phi,eps,totNFM,dim,&
                  ng,Nmat,order,nregion,it,parametre)
        implicit none
        integer(kind=4), intent(in) :: ng,dim,totNFM,Nmat,order,nregion,it
        real(kind=8), dimension(Nmat,ng), intent(in) :: SigT,NusigF,Chi
        real(kind=8), dimension(Nmat,ng,ng), intent(in) :: SigS
        real(kind=8), dimension(dim), intent(in) :: phi
        real(kind=8), dimension(nregion), intent(in) :: dcell
        real(kind=8), dimension(2), intent(in) :: albedo
        CHARACTER(50), intent(in) :: start,tm
        real(kind=8), intent(in) :: eps,k_eff
        integer(kind=4), intent(in) :: parametre
        integer(kind=4) :: i,j
        open (100,file='app/Output/OUTPUT_CP.TXT')
        write (100, FMT=* ) '********************************************************************************'
        write (100, FMT=* ) 'ERSN, UNIVERSITY ABDELMALEK ESSAADI FACULTY OF SCIENCES - TETOUAN, MOROCCO'
        write (100, FMT=* ) 'CODE  DEVELOPED  BY  MOHAMED  LAHDOUR,  PHD  STUDENT'
        write (100, FMT=* ) 'OpenRSN:         CP COLLISION PROBABILITY METHOD'
        write (100, FMT=* ) 'DIMENSION:       ONE DIMENSION (1D) '
        if     (parametre == 1) then 
        write (100, FMT=* ) 'GEOMETRY:        SLAB'
        elseif (parametre == 2) then 
        write (100, FMT=* ) 'GEOMETRY:        CYLINDER'
        elseif (parametre == 3) then 
        write (100, FMT=* ) 'GEOMETRY:        SPHERE'
        endif
        write (100, FMT=* ) 'VERSION NUMBER:  1.2'
        write (100, FMT=* ) 'VERSION DATE:    20  August  2019'
        write (100,3010) 'RAN ON:          ', start,'(H:M:S)'
        write (100, FMT=* ) '********************************************************************************'
        write (100, FMT=* ) '           ----------------------------------------------------------' 
        write (100, FMT=* ) '                     INPUT  PARAMETER - VALUES  FROM  INPUT'              
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'ENERGY GROUPS NUMBER:                    ',ng
        write (100, FMT=* ) 'REGIONS NUMBER:                          ',nregion
        write (100, FMT=* ) 'MATERIALS NUMBER:                        ',Nmat
        if     (parametre == 1) then 
        write (100,3040)    'SIZE FOR EACH REGION PER [CM]:           ',dcell  
        else
        write (100,3040)    'RAY FOR EACH REGION PER [CM]:            ',dcell
        endif     
        write (100, FMT=* ) 'L-ORDER LEGENDRE POLYNOMIAL:             ',order-1
        write (100, FMT=* ) 'TOTAL NUMBER OF FINE MESHES:             ',totNFM
        write (100,3050)    'CONVERGENCE CRITERION of KEFF AND FLUX:  ',eps
        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '                      CALCULATION  RUN-TIME  PARAMETERS  CP' 
        write (100, FMT=* ) '           ----------------------------------------------------------'
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
            write(100,3080) j,SigT(i,j),SigT(i,j)-SigS(i,j,j),NusigF(i,j),SigS(i,j,j),Chi(i,j)
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
        if  ( albedo(1) == 1.0D0 .and. albedo(2) == 1.0D0 ) then
        write (100,3090)    'K-INF                    =',k_eff
        else
        write (100,3090)    'K-EFF                    =',k_eff
        endif
        write (100,3020)    'N. OUTER ITERATIONS      =',it
        write (100,4000)    'TOTAL EXECUTION TIME     =',tm,'(H:M:S)'
        write (100, FMT=* ) ''
        write (100, FMT=* ) '********************************************************************************'
        2000 format(1x,1p,i11,5x,200e16.5) 
        3000 format(1x,A14,2x,300(A14,i2))  
        3010 format(1x,A17,A22,A10)
        3020 format(1x,A26,4x,i10)
        3040 format(1x,A33,2x,200F10.5)
        3050 format(1x,1p,A41,4x,e8.1)
        3070 format(1x,A18,i4)
        3080 format(1x,1p,i11,5x,e16.5,e16.5,e16.5,e16.5,e16.5)
        3090 format(1x,A26,6x,f8.6)
        4000 format(1x,A26,4x,A10,A10)
        close(100)
end subroutine Output
!BL26
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
!BL27
    subroutine title1(parametre) 
       integer(kind=4), intent(in) :: parametre  
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
       print*, 'Method   | The collision probability method (CP)            '  
       print*, 'Dimension| One dimensions (1D)                              '
       if     (parametre == 1) then 
       print*, 'Geometry | Slab                                             '  
       elseif (parametre == 2) then 
       print*, 'Geometry | Cylinder                                         ' 
       elseif (parametre == 3) then 
       print*, 'Geometry | Sphere                                           ' 
       endif
       print*, '____________________________________________________________'
       print*, ''
    end subroutine title1
!BL23
    subroutine title2()
       write(*,FMT=*)'**********************************************************'
       write(*,FMT=*)'                         Finished                         '                             
       write(*,FMT=*)'**********************************************************'  
    end subroutine title2


 
