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
end subroutine inverse
! ------------------------------------------------------------------------------
 subroutine invers(c,b,n)
 ! subroutine to calculate the inverse of a matrix using Gauss-Jordan elimination
 ! the inverse of matrix a(n,n) is calculated and stored in the matrix b(n,n)
INTEGER         ,INTENT(in)    :: n
DOUBLE PRECISION,INTENT(in)  ::  c(n,n)
DOUBLE PRECISION,INTENT(out) ::  b(n,n)
integer :: i,j,k,l,m,irow
DOUBLE PRECISION :: big,dum,a(n,n)
!build the identity matrix
a = c
 do i = 1,n
 do j = 1,n
 b(i,j) = 0.0
 end do
 b(i,i) = 1.0
 end do

 do i = 1,n ! this is the big loop over all the columns of a(n,n)
 ! in case the entry a(i,i) is zero, we need to find a good pivot; this pivot
 ! is chosen as the largest value on the column i from a(j,i) with j = 1,n
 big = a(i,i)
 do j = i,n
 if (a(j,i).gt.big) then
 big = a(j,i)
 irow = j
 end if
 end do
 ! interchange lines i with irow for both a() and b() matrices
 if (big.gt.a(i,i)) then
 do k = 1,n
 dum = a(i,k) ! matrix a()
 a(i,k) = a(irow,k)
 a(irow,k) = dum
 dum = b(i,k) ! matrix b()
 b(i,k) = b(irow,k)
 b(irow,k) = dum
 end do
 end if
 ! divide all entries in line i from a(i,j) by the value a(i,i);
 ! same operation for the identity matrix
 dum = a(i,i)
 do j = 1,n
 a(i,j) = a(i,j)/dum
 b(i,j) = b(i,j)/dum
 end do
 ! make zero all entries in the column a(j,i); same operation for indent()
 do j = i+1,n
 dum = a(j,i)
 do k = 1,n
 a(j,k) = a(j,k) - dum*a(i,k)
 b(j,k) = b(j,k) - dum*b(i,k)
 end do
 end do
 end do

 ! substract appropiate multiple of row j from row j-1
 do i = 1,n-1
 do j = i+1,n
 dum = a(i,j)
 do l = 1,n
 a(i,l) = a(i,l)-dum*a(j,l)
 b(i,l) = b(i,l)-dum*b(j,l)
 end do
 end do
 end do

 end subroutine invers
!BL
    subroutine Aitken_SOR(A,phi_li,F,phi,eps,ng,totNFM,Norder,dim,ddim)
!      Coupling of Aitken's delta-squared with Successive Over-Relaxation  
!      -----------------------------------------------------------------
       implicit none
       ! global variables
       integer(kind=4), intent(in) :: ng,totNFM,dim,ddim,Norder
       real(kind=8), dimension(ddim,ddim), intent(in) :: A
       real(kind=8), dimension(ddim,ddim), intent(in) :: F
       real(kind=8), dimension(ddim), intent(in) :: phi_li
       real(kind=8), intent(in) :: eps
       real(kind=8), dimension(dim), intent(out) :: phi
       ! locall variables
       real(kind=8), dimension(ddim,ddim) :: ai,L,U,AA,dinv,INI,ai1,ai2,ainv
       real(kind=8), dimension(ddim) :: x0,x1,x2,phi_li0,phi_li1
       real(kind=8), dimension(ddim) :: sold
       real(kind=8) :: err1,err2,sumnew,sumold
       real(kind=8) ::  err_keff,err_phi, k_eff=1.0,k_eff0=1.0,norme
       integer(kind=4) ::  iter = 0,i,j,k1,k0
       real(kind=8), dimension(ng) :: moy 
       err_keff = 1.0
       ai=0.;L=0.;U=0.;AA=0.;dinv=0.;INI=0.;ai1=0.;ai2=0.
       !--------------------------------
       do i = 1,ddim
          dinv(i,i) = 1./A(i,i)
          INI(i,i) = 1.0
       enddo
       AA = matmul(dinv,A)
       do i = 1,ddim
          do j = i,ddim
             if (i.ne.j) then
             L(j,i) = -AA(j,i)
             U(i,j) = -AA(i,j)
             endif
          enddo
       enddo
       ai1 = INI - L
       call  inverse(ai1,ainv,ddim) !matinv(ai1,ainv,ddim)  invers(ai1,ainv,ddim) 
       ai  = matmul(dinv,F)
       ai2 = matmul(ainv,ai)
       ai1 = matmul(ainv,U)
       !--------------------------------
       x0 = phi_li
       do i = 1,2
           sold = matmul(ai2,x0)/k_eff
           sumold = sum(matmul(F,x0))
           if (i==1) then
           x1   =  k_eff*matmul(ai1,x0) + sold 
           x0 = x1
           else
           x2   =  k_eff*matmul(ai1,x0) + sold
           x0 = x2
           endif
           sumnew = sum(matmul(F,x0))
           k_eff0 = k_eff
           k_eff =  k_eff*(sumnew/sumold)
       enddo
       !--------------------------------
           phi_li0 = k_eff*x2 - k_eff0*x1
       ! Outer Iteration
       do while ( err_keff >= eps )
           err_phi = 1.0
           iter = iter + 1
           sold = matmul(ai2,phi_li0)/k_eff
           sumold = sum(matmul(F,phi_li0))
           if ( iter > 1000 ) then
           print*, 'error(unable to converge(1))'
           stop
           end if
               ! Inner Iteration
               do while ( err_phi>= eps )
               phi_li1 = matmul(ai1,phi_li0) + sold 
               err1   = maxval(abs(phi_li1))
               err2   = maxval(abs(phi_li1-phi_li0))
!              critere de convergence sur les vecteurs porpres   
               err_phi = err2/err1 
               phi_li0 =  phi_li1
               enddo
           sumnew = sum(matmul(F,phi_li1))
!          critere de convergence sur les valeurs porpres
           k_eff =  k_eff*(sumnew/sumold)
           err1  = abs(k_eff-k_eff0)
           err_keff = err1/k_eff
           norme = sqrt(sum(phi_li1*phi_li1))
           phi_li0 = phi_li0/norme
           x1 = x2
           x2 = phi_li0     
           write(*,2000)iter,k_eff,err_keff
           k_eff0 = k_eff
       enddo
           phi_li1 = x1 - phi_li0/(1.+ k_eff)
      !Normalized flux --------------------------
             moy = 0.0
             k1 = 1
             j = 1
             do k0 = 1,ng
                 do i=k1,(k1-1)+totNFM
                    moy(k0) = moy(k0) + phi_li1(i)
                 end do
                 do i=k1,(k1-1)+totNFM
                    phi(j) = phi_li1(i)/(moy(k0)/totNFM)
                    j = j + 1
                 end do
                 k1 = k1 + totNFM*(Norder+1)/2
             enddo
       ! Format
       2000 format(t3,"Iteration",i4,":",5x,"===>",5x,"keff =",F9.6,5x,&
                  "===>",5x,"res =",e10.3)
    end subroutine Aitken_SOR
!BL
    subroutine GSL_SOR(A,phi_li,F,ML,MU,phi,eps,ng,totNFM,Norder,dim,ddim)
       implicit none
       ! global variables
       integer(kind=4), intent(in) :: ng,totNFM,dim,ddim,Norder
       real(kind=8), dimension(ddim,ddim), intent(in) :: A
       real(kind=8), dimension(ddim,ddim), intent(in) :: F,ML,MU
       real(kind=8), dimension(ddim), intent(in) :: phi_li
       real(kind=8), intent(in) :: eps
       real(kind=8), dimension(dim), intent(out) :: phi
       ! locall variables
       real(kind=8), dimension(ddim,ddim) :: U,L,ai,ainv,dinv,AA,INI,ai1,ai2,B
       real(kind=8), dimension(ddim) :: phi_li0,phi_li1,sold
       real(kind=8), dimension(ng) :: moy
       real(kind=8) :: sumold,sumnew,norme,omega=1.4
       real(kind=8) ::  k_eff=1.,k_eff0=1.,err_keff,err_phi,err1,err2
       integer(kind=4) :: i,j,iter=0,iner=0,k0,k1
       INI = 0.0;U=0.;L=0.;AA=0.;ai1=0.;ai=0.;ai2=0.;norme=0.0
       err_keff = 1.0
       !----------------------
       do i = 1,ddim
          dinv(i,i) = 1./A(i,i)
          INI(i,i) = 1.0
       enddo
       B = A - ML - MU 
       AA = matmul(dinv,B) 
       ai = matmul(dinv,F)
       do i = 1,ddim
          do j = i,ddim
             if (i.ne.j) then
             L(j,i) = -AA(j,i)
             U(i,j) = -AA(i,j)
             endif
          enddo
       enddo
       ai2 = INI - omega*L
       call inverse(ai2,ainv,ddim)
       ai  = omega*matmul(ainv,ai)
       ai1 = matmul(ainv,omega*U+(1.-omega)*INI)
       phi_li0 = phi_li
       ! Outer Iteration
       do while ( err_keff >= eps )
                err_phi =1.
                iter = iter + 1
                sold = matmul(ai,phi_li0)/k_eff
                sumold = sum(matmul(F,phi_li0))
                if ( iter > 1000 ) then
                print*, 'error(unable to converge(1))'
                stop
                end if
               ! Inner Iteration
               do while ( err_phi>= eps )
                         phi_li1 = matmul(ai1,phi_li0) + sold 
!                        critere de convergence sur les vecteurs porpres
                         err1   = maxval(abs(phi_li1))
                         err2   = maxval(abs(phi_li1-phi_li0))
                         err_phi = err2/err1 
                         phi_li0 = phi_li1
               enddo
                sumnew = sum(matmul(F,phi_li1))  
                k_eff =  k_eff*(sumnew/sumold)
                err1 = abs(k_eff-k_eff0)
                err_keff = err1/k_eff
                k_eff0 = k_eff 
                norme = sqrt(sum(phi_li1*phi_li1))
                phi_li0 = phi_li1/norme 
                write(*,2000)iter,k_eff,err_keff                       
       enddo
       !Normalized flux --------------------------
             moy = 0.0
             k1 = 1
             j = 1
             do k0 = 1,ng
                 do i=k1,(k1-1)+totNFM
                    moy(k0) = moy(k0) + phi_li1(i)
                 end do
                 do i=k1,(k1-1)+totNFM
                    phi(j) = phi_li1(i)/(moy(k0)/totNFM)
                    j = j + 1
                 end do
                 k1 = k1 + totNFM*(Norder+1)/2
             enddo
       ! Format
       2000 format(t3,"Iteration",i4,":",5x,"===>",5x,"keff =",F9.6,5x,&
                  "===>",5x,"res =",e10.3)
    end subroutine GSL_SOR
!BL 
   subroutine Aitken(ainv,A,phi_li,F,Delta,phi,eps,ng,totNFM,Norder,dim,ddim)
!      -------------------  The inverse power method  -------------------
!                           (A - 1/Keff*F)*phi = 0
! 
!      -----------------------------------------------------------------
       implicit none
       ! global variables
       integer(kind=4), intent(in) :: ng,totNFM,dim,ddim,Norder
       real(kind=8), dimension(ddim,ddim), intent(in) :: ainv,A
       real(kind=8), dimension(ddim,ddim), intent(in) :: F
       real(kind=8), dimension(ddim), intent(in) :: phi_li
       real(kind=8), dimension(totNFM), intent(in) :: Delta
       real(kind=8), intent(in) :: eps
       real(kind=8), dimension(dim), intent(out) :: phi
       ! locall variables
       real(kind=8), dimension(ddim,ddim) :: ai
       real(kind=8), dimension(ddim) :: x0,x1,x2,phi_li0,phi_li1
       real(kind=8), dimension(ddim) :: gar,vec1,vec2
       real(kind=8) :: s1,s2,err1,err2
       real(kind=8) ::  err_k_eff,err_phi,eval1=0.0, eval2,test
       integer(kind=4) ::  iter = 0,i,j,k1,k0
       real(kind=8), dimension(ng) :: moy 
       err_k_eff = 1.0
       err_phi = 1.0
       x0 = phi_li
       ai = matmul(ainv,F)
       do i = 1,2
           gar    = matmul(ai,x0)
           vec1   = matmul(A,x0)
           vec2   = matmul(F,x0)
           s1     = dot_product(vec1,vec2)
           s2     = dot_product(vec2,vec2)
           eval2  = s1/s2
           eval1  = eval2
           if (i==1) then
           x1   = gar*eval1 
           x0 = x1
           else
           x2   = gar*eval1 
           endif
       enddo
           x0 = phi_li
           phi_li0 = eval2*x2 - eval1*x1
           eval1 = 0.
       ! Outer Iteration
       do while ( err_k_eff>=eps .and. err_phi>= eps )
           iter = iter + 1
           if ( iter > 1000 ) then
           print*, 'error(unable to converge(1))'
           stop
           end if
           !phi_li0 = phi_li0/(sqrt(sum(phi_li0*phi_li0)))
           !gar    = matmul(ai,phi_li0)
           vec1   = matmul(A,phi_li0)
           vec2   = matmul(F,phi_li0)
           s1     = dot_product(vec1,vec2)
           s2     = dot_product(vec2,vec2)
           eval2  = s1/s2
!          critere de convergence sur les valeurs porpres
           err_k_eff =  abs(eval2-eval1)/eval2
           eval1  = eval2
           x0 = x1; x1 = x2; x2 = phi_li0
           phi_li1   =  x0 - (x1-x0)*(x1-x0)/(x2-2*x1+x0)
           !phi_li1 = phi_li1/(sqrt(sum(phi_li1*phi_li1)))
           !phi_li1 = eval1*gar
           err1   = maxval(abs(phi_li1))
           err2   = maxval(abs(phi_li1-phi_li0))
!          critere de convergence sur les vecteurs porpres   
           err_phi = err2/err1       
           write(*,2000)iter,1/eval1,err_k_eff
           phi_li0   = phi_li1
       enddo
       !Normalized flux --------------------------
             moy = 0.0
             k1 = 1
             j = 1
             do k0 = 1,ng
                 do i=k1,(k1-1)+totNFM
                    moy(k0) = moy(k0) + phi_li1(i)
                 end do
                 do i=k1,(k1-1)+totNFM
                    phi(j) = phi_li1(i)/(moy(k0)/totNFM)
                    j = j + 1
                 end do
                 k1 = k1 + totNFM*(Norder+1)/2
             enddo
       ! Format
       2000 format(t3,"Iteration",i4,":",5x,"===>",5x,"keff =",F9.6,5x,&
                  "===>",5x,"res =",e10.3)
    end subroutine Aitken
!BL
subroutine Decomp_QR(n,A,Q,R)
integer(kind=4), intent(in) :: n
real(kind=8), dimension(n,n), intent(inout) :: A
real(kind=8), dimension(n,n), intent(inout) :: Q,R
real(kind=8), dimension(n,n) :: INI,H
real(kind=8), dimension(n) :: v
real(kind=8) :: alpha,betha,c
integer(kind=4) :: i,j,k
INI=0.;v=0.;H=0
   do k0 =1,n
      INI(k0,k0) = 1.
   enddo
H=INI
do k = 1,n-1
   alpha = sqrt(sum(A(k:n,k)*A(k:n,k)))
   betha = alpha**2 - alpha*A(k,k)
   v(k) = A(k,k) - alpha
   do i = k+1,n
      v(i) = A(i,k)
   enddo
   do j = k,n
      c = 1./betha*(sum(v(k:n)*A(k:n,j)))
      do i = k,n
         A(i,j) = A(i,j) - c*v(i)
      enddo
   enddo
   do j=1,n
      c = 1/betha*(sum(v(k:n)*H(k:n,j)))
      do i = k,n
         H(i,j) = H(i,j) - c*v(i)
      enddo
   enddo
enddo
Q = transpose(H)
R = A
end subroutine Decomp_QR   

!BL
       subroutine GSL(A,phi_li,F,phi,eps,ng,totNFM,Norder,dim,ddim)
!      -------------------  The inverse power method  ------------------
!                           (A - 1/Keff*F)*phi = 0
! 
!      -----------------------------------------------------------------
       implicit none
       ! global variables
       integer(kind=4), intent(in) :: ng,totNFM,dim,ddim,Norder
       real(kind=8), dimension(ddim,ddim), intent(in) :: A
       real(kind=8), dimension(ddim,ddim), intent(in) :: F
       real(kind=8), dimension(ddim), intent(in) :: phi_li
       real(kind=8), intent(in) :: eps
       real(kind=8), dimension(dim), intent(out) :: phi
       ! locall variables
       real(kind=8), dimension(ddim,ddim) :: U,L,ai,ainv,dinv,AA,INI,ai1,ai2
       real(kind=8), dimension(ddim) :: phi_li0,phi_li1,sold
       real(kind=8), dimension(ng) :: moy
       real(kind=8) :: sumold,sumnew,norme
       real(kind=8) ::  k_eff=1.,k_eff0=1.,err_keff,err_phi,err1,err2
       integer(kind=4) :: i,j,iter=0,iner=0,k0,k1
       INI = 0.0;U=0.;L=0.;AA=0.;ai1=0.;ai=0.;ai2=0.;norme=0.0
       err_keff = 1.0
       !----------------------
       do i = 1,ddim
          dinv(i,i) = 1./A(i,i)
          INI(i,i) = 1.0
       enddo
       AA = matmul(dinv,A)
       ai = matmul(dinv,F)
       do i = 1,ddim
          do j = i,ddim
             if (i.ne.j) then
             L(j,i) = -AA(j,i)
             U(i,j) = -AA(i,j)
             endif
          enddo
       enddo
       ai2 = INI - L
       call inverse(ai2,ainv,ddim)
       ai  = matmul(ainv,ai)
       ai1 = matmul(ainv,U)
       norme = sqrt(sum(phi_li(:)**2))
       phi_li0 = phi_li/norme
       phi_li1 = phi_li0
       ! Outer Iteration
       do iter = 1,1000
                sold = matmul(ai,phi_li1)/k_eff
                sumold = sum(matmul(F,phi_li1))
                do iner = 1,1000
                         phi_li1 = matmul(ai1,phi_li1) + sold  
!             critere de convergence sur les vecteurs porpres
                         err1   = maxval(abs(phi_li1))
                         err2   = maxval(abs(phi_li1-phi_li0))
                         err_phi = err2/err1 
                         phi_li0 = phi_li1
                        ! print*,iner,err_phi
                         if  (err_phi <= eps) exit
                enddo
                sumnew = sum(matmul(F,phi_li1))  
                k_eff =  k_eff*(sumnew/sumold)
                err1 = abs(k_eff-k_eff0)
                err_keff = err1/k_eff
                k_eff0 = k_eff 
                norme = sqrt(sum(phi_li1(:)**2))
                phi_li1 = phi_li1/norme 
                write(*,2000)iter,k_eff,err_keff
                if  (err_keff <= eps) exit                       
       enddo
       !------------------------------------------
       !Normalized flux --------------------------
             moy = 0.0
             k1 = 1
             j = 1
             do k0 = 1,ng
                 do i=k1,(k1-1)+totNFM
                    moy(k0) = moy(k0) + phi_li1(i)
                 end do
                 do i=k1,(k1-1)+totNFM
                    phi(j) = phi_li1(i)/(moy(k0)/totNFM)
                    j = j + 1
                 end do
                 k1 = k1 + totNFM*(Norder+1)/2
             enddo
       ! Format
       2000 format(t3,"Iteration",i4,":",5x,"===>",5x,"keff =",F9.6,5x,&
                  "===>",5x,"res =",e10.3)
       end subroutine GSL
!BL
       subroutine Solve_QR(A,F,eps,phi,phi_li,ng,totNFM,Norder,ddim)
!      -------------------  Solve QR Decomposition  -------------------
!                           (A - 1/Keff*B)*phi = 0
! 
!      ------------------------------------------------------------------
       implicit none
       integer(kind=4), intent(in) :: ng,totNFM,ddim,Norder
       real(kind=8), intent(in) :: eps
       real(kind=8), dimension(ddim,ddim), intent(in) :: A,F
       real(kind=8), dimension(ddim), intent(in) :: phi_li
       real(kind=8), dimension(ng*totNFM), intent(out) :: phi
       real(kind=8), dimension(ddim,ddim) :: D,INI,U,dinv,AA,ai1,ai2,Q,R,W
       real(kind=8), dimension(ddim) :: sold
       real(kind=8), dimension(ddim) :: phi1,phi2
       real(kind=8) :: sumold,sumnew,norme
       real(kind=8), dimension(ng) :: moy
       real(kind=8) :: err1,err2,test,err_keff=1.0,err_phi,k_eff=1.,k_eff0=1.
       integer(kind=4) :: i,j,k0,k1,iner,iter
       D=0.;U=0.;INI=0.;dinv=0.;AA=0.;ai1=0.;ai2=0.
       W = A
       call  Decomp_QR(ddim,W,Q,R)
       do i = 1,ddim
          INI(i,i) = 1.
          D(i,i) = R(i,i)
       enddo
       call matinv(D,dinv,ddim)
       AA = matmul(dinv,R)
       U = INI-AA
       iter = 0
       !phi1 = phi_li
       norme = sqrt(sum(phi_li*phi_li))
       phi1 = phi_li/norme 
       print*,phi1
       ai1 =  matmul(transpose(Q),F) 
       ai2 =  matmul(dinv,ai1)
       do while ( err_keff >= eps) 
           iner = 0
           err_phi = 1.0      
           iter = iter + 1
           if ( iter > 1000 ) then
           print*, 'error(unable to converge(1))'
           stop
           end if
           sold   = matmul(ai2,phi1)/k_eff
           sumold = sum(matmul(F,phi1))
           do while ( err_phi >= eps )
                iner = iner + 1
                phi2  = matmul(U,phi1) + sold 
!               critere de convergence sur les valeurs porpres                       
                err1   = maxval(abs(phi2))
                err2   = maxval(abs(phi2-phi1))
!               critere de convergence sur les vecteurs porpres   
                err_phi = err2/err1       
                phi1   = phi2 
                print*,iner,err_phi
           enddo
           print*,phi1
           sumnew = sum(matmul(F,phi1))  
           k_eff =  k_eff*(sumnew/sumold)
           err1 = abs(k_eff-k_eff0)
           err_keff = err1/k_eff
           k_eff0 = k_eff 
           norme = sqrt(sum(phi1*phi1))
           phi1 = phi1/norme    
           if (iter == 1)  then
                   test = err_keff
           elseif (iter >= 10 .and. err_keff > test) then
           print*, 'error(unable to converge(2) )'
           print*, 'because error in iteration ten is sup to one'
           stop
           end if
           write(*,2000)iter,k_eff,err_keff
       end do
       2000 format(t3,"Iteration",i4,":",5x,"===>",5x,"keff =",F9.6,5x,"===>",5x,"res =",e10.3)
       !Normalized flux --------------------------
             moy = 0.0
             k1 = 1
             j = 1
             do k0 = 1,ng
                 do i=k1,(k1-1)+totNFM
                    moy(k0) = moy(k0) + phi1(i)
                 end do
                 do i=k1,(k1-1)+totNFM
                    phi(j) = phi1(i)/(moy(k0)/totNFM)
                    j = j + 1
                 end do
                 k1 = k1 + totNFM*(Norder+1)/2
             enddo

       end subroutine Solve_QR
!BL
subroutine conjgrad(ng,np,totNFM,A,F,Delta,x,phi,eps)  ! Solve Ax=b with OrthoMin algo
!  use constants
  implicit none
  integer, intent(in)                :: ng,np,totNFM
  real, dimension(np,np), intent(in) :: A, F
  real, dimension(totNFM), intent(in)    :: Delta
  real, dimension(np), intent(in) :: x
  real, dimension(totNFM), intent(out) :: phi
  real , intent(in) :: eps
  real, dimension(np,np)  :: B
  real, dimension(np) :: xi,vi,ri,s, vec1, vec2
  real, dimension(totNFM) :: moy
  real                :: di=1.,dip,wi,psi,d0,eval1=1.0, eval2,s1,s2,err_k_eff=1.0
  integer i,j,k0,k1
  B = 0
  ! ------------------------
  do i = 1,totNFM
     B(i,i) = Delta(i)*F(i,i)
  enddo
  ! ------------------------
  ! ... initial state ...
  j = 0
  xi  = x  ! initial guess
  ! ... CG MINRES iteration ...
      do while ( err_k_eff >= eps )
     !.........
      s = eval1*matmul(B,xi)
      ri  = s - MATMUL(A,xi)
      vi  = ri
      di  = DOT_PRODUCT(ri,ri)
      dip = DOT_PRODUCT(vi,MATMUL(A,vi))
     !.........
           do while ( di >= 1.e-6)
           wi  = di/dip
           ri  = ri - MATMUL(A,vi)*wi
           xi  = xi + vi*wi
           d0  = di
           di  = DOT_PRODUCT(ri,ri)
           psi = -di/d0
           vi  = ri-vi*psi
           dip = DOT_PRODUCT(vi,MATMUL(A,vi))
           enddo
      vec1   = matmul(A,xi)
      vec2   = matmul(B,xi)
      s1     = dot_product(vec1,vec2)
      s2     = dot_product(vec2,vec2)
      eval2  = s1/s2
      err_k_eff =  abs(eval2-eval1)/eval2
      eval1  = eval2
      s1 = s2 
      j   = j + 1
      write(*,2000)j,1/eval1,err_k_eff
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
!BL1
    subroutine matinv(a,ainv,n)
       implicit none
       integer(kind=4), intent(in) :: n
       real(kind=8), dimension(n,n), intent(in) :: a
       real(kind=8), dimension(n,n), intent(out) :: ainv
       real(kind=8), dimension(n,2*n) :: b
       integer(kind=8) :: i,j,k
       real(kind=8) :: pivot=0.0,xnum
       open (100,file='app/Output/test.h')
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
       !write(*,*)'Invers Matrx A'
       !write(100,'(200f8.4)') transpose(ainv)
       return
    end subroutine matinv 
!BL2
    subroutine aleig2(A,phi_li,F,phi,eps,ng,totNFM,dim,ddim)
!      -------------------  The inverse power method  ------------------
!                           (A - 1/Keff*F)*phi = 0
! 
!      -----------------------------------------------------------------
       implicit none
       ! global variables
       integer(kind=4), intent(in) :: ng,totNFM,dim,ddim
       real(kind=8), dimension(ddim,ddim), intent(in) :: A
       real(kind=8), dimension(ddim,ddim), intent(in) :: F
       real(kind=8), dimension(ddim), intent(in) :: phi_li
       real(kind=8), intent(in) :: eps
       real(kind=8), dimension(dim), intent(out) :: phi
       ! locall variables
       real(kind=8), dimension(ddim,ddim) :: D,O
       real(kind=8), dimension(ddim) :: phi_li0,phi_li1,sold,snew
       real(kind=8) :: sumold,sumnew  !,norme
       real(kind=8) ::  k_eff=1.,k_eff0=1.,norme,err_keff,err_phi,err1,err2
       integer(kind=4) :: i,k1,k0,iter=0,iner=0
       real(kind=8), dimension(ng) :: moy 
       err_keff = 1.0
       phi_li0 = phi_li
       !----------------------
       D=0.; O= 0.
       do i = 1,ddim
          D(i,i) = A(i,i)
       enddo 
       O = A - D 
       ! Outer Iteration
       do iter = 1,1
                sold = matmul(F,phi_li0)
                sumold= sum(sold)
                sold = sold/(sumold/ddim)
                sumold = sum(sold)
                do iner = 1,100 
                         phi_li1 = sold/k_eff - matmul(O,phi_li0)
                         do i = 1,ddim
                            phi_li1(i) = phi_li1(i)/D(i,i)
                         enddo  
                         err1   = maxval(abs(phi_li1))
                         err2   = maxval(abs(phi_li1- phi_li0))
!             critere de convergence sur les vecteurs porpres   
                         err_phi = err2/err1  
                         norme = sqrt(sum(phi_li1(:)*phi_li1(:)))
                         phi_li0 = phi_li1/norme
                         print*,iner,sum(phi_li1)
                         if  (err_phi <= eps) exit
                enddo
                snew = matmul(F,phi_li1) 
                print*,snew
                sumnew = sum(snew)  
                snew = snew/(sumnew/totNFM) 
                sumnew = sum(snew)
                k_eff = k_eff*sumnew/sumold
                err_keff = abs(k_eff-k_eff0)/k_eff
                k_eff0=k_eff  
                if  (err_keff <= eps) exit                       
           write(*,2000)iter,k_eff,err_keff
       enddo
       !------------------------------------------
       !Normalized flux --------------------------
             moy = 0.0
             k1 = 1
             do k0 = 1,ng
                 do i=k1,totNFM*k0
                    moy(k0) = moy(k0) + phi_li1(i)
                 end do
                 do i=k1,totNFM*k0
                    phi(i) = phi_li1(i)/(moy(k0)/totNFM)
                 end do
                 k1 = k1 + totNFM
             enddo
       ! Format
       2000 format(t3,"Iteration",i4,":",5x,"===>",5x,"keff =",F9.6,5x,&
                  "===>",5x,"res =",e10.3)
    end subroutine aleig2
!BL1
    subroutine LU_Factorisation(L,U,phi_li,F,phi,eps,ng,totNFM,Norder,iter,inter,k_eff,dim,ddim)
!      -------------------  The inverse power method  ------------------
!                           (A - 1/Keff*F)*phi = 0
! 
!      -----------------------------------------------------------------
       implicit none
       ! global variables
       integer(kind=4), intent(in) :: ng,totNFM,dim,ddim,Norder
       real(kind=8), dimension(ddim,ddim), intent(in) :: L,U
       real(kind=8), dimension(ddim,ddim), intent(in) :: F
       real(kind=8), dimension(ddim), intent(in) :: phi_li
       real(kind=8), intent(in) :: eps
       integer(kind=4), intent(out) :: iter,inter
       real(kind=8), intent(out) :: k_eff
       real(kind=8), dimension(dim), intent(out) :: phi
       ! locall variables
       real(kind=8), dimension(ddim) :: phi_li0,phi_li1,Y,sold,snew
       real(kind=8) :: err1,err2,sumnew,sumold!,norme
       real(kind=8) ::  err_k_eff,err_phi,k_eff0=1.
       integer(kind=4) ::  i,j,k1,k0
       real(kind=8), dimension(ng) :: moy 
       iter = 0
       inter = 0
       err_k_eff = 1.0
       phi_li0 = phi_li 
       sold = matmul(F,phi_li0) 
       sumold=sum(sold)
       sold = sold/(sumold/totNFM)
       sumold = sum(sold)
       ! Outer Iteration
       do while ( err_k_eff>= eps)
           err_phi = 1.0
           iter = iter + 1
           call ForSub(ddim,L,sold,Y)
           if ( iter > 1000 ) then
           print*, 'error(unable to converge(1))'
           stop
           end if 
           do while ( err_phi >= eps) 
              call BackSub(ddim,U,Y,phi_li1)
              err1   = maxval(abs(phi_li1))
              err2   = maxval(abs(phi_li1- phi_li0))
!             critere de convergence sur les vecteurs porpres   
              err_phi = err2/err1  
              phi_li0   = phi_li1
              inter = inter + 1 
           enddo 
           !norme = sqrt(sum(phi_li1*phi_li1))
           !phi_li1 = phi_li1/norme     
           snew = matmul(F,phi_li1)
           sumnew=sum(snew)
           k_eff = sumnew/sumold
           sold = snew/(sumnew/totNFM)
           sumold = sum(sold)
           err_k_eff = abs(k_eff-k_eff0)/k_eff
!          critere de convergence sur les valeurs porpres                             
           write(*,2000)iter,k_eff,err_k_eff
           k_eff0 = k_eff
       enddo
       !------------------------------------------
       !Normalized flux --------------------------
             moy = 0.0
             k1 = 1
             j = 1
             do k0 = 1,ng
                 do i=k1,(k1-1)+totNFM
                    moy(k0) = moy(k0) + phi_li1(i)
                 end do
                 do i=k1,(k1-1)+totNFM
                    phi(j) = phi_li1(i)/(moy(k0)/totNFM)
                    j = j + 1
                 end do
                 k1 = k1 + totNFM*(Norder+1)/2
             enddo
       ! Format
       2000 format(t3,"Iteration",i4,":",5x,"===>",5x,"keff =",F9.6,5x,&
                  "===>",5x,"res =",e10.3)
    end subroutine LU_Factorisation
!BL2
    subroutine aleig(ainv,A,phi_li,F,phi,eps,ng,totNFM,Norder,dim,ddim)
!      -------------------  The inverse power method  -------------------
!                           (A - 1/Keff*F)*phi = 0
! 
!      -----------------------------------------------------------------
       implicit none
       ! global variables
       integer(kind=4), intent(in) :: ng,totNFM,dim,ddim,Norder
       real(kind=8), dimension(ddim,ddim), intent(in) :: ainv,A
       real(kind=8), dimension(ddim,ddim), intent(in) :: F
       real(kind=8), dimension(ddim), intent(in) :: phi_li
       real(kind=8), intent(in) :: eps
       real(kind=8), dimension(dim), intent(out) :: phi
       ! locall variables
       real(kind=8), dimension(ddim,ddim) :: ai
       real(kind=8), dimension(ddim) :: phi_li0,phi_li1
       real(kind=8), dimension(ddim) :: gar,vec1,vec2
       real(kind=8) :: s1,s2,err1,err2
       real(kind=8) ::  err_k_eff,err_phi,eval1=0.0, eval2
       integer(kind=4) ::  iter = 0,i,k1,k0,j
       real(kind=8), dimension(ng) :: moy 
       err_k_eff = 1.0
       err_phi = 1.0
       phi_li0 = phi_li
       ai = matmul(ainv,F)  
       ! Outer Iteration
       do while ( err_k_eff>= eps .and. err_phi >= eps)
           iter = iter + 1
           if ( iter > 1000 ) then
           print*, 'error(unable to converge(1))'
           stop
           end if
           gar    = matmul(ai,phi_li0)
           vec1   = matmul(A,phi_li0)
           vec2   = matmul(F,phi_li0)
           s1     = dot_product(vec1,vec2)
           s2     = dot_product(vec2,vec2)
           eval2  = s1/s2
!          critere de convergence sur les valeurs porpres               
           err_k_eff =  abs(eval2-eval1)/eval2
           eval1  = eval2 
           phi_li1   = gar*eval1         
           err1   = maxval(abs(phi_li1))
           err2   = maxval(abs(phi_li1-phi_li0))
!          critere de convergence sur les vecteurs porpres   
           err_phi = err2/err1       
           phi_li0   = phi_li1 
           write(*,2000)iter,1/eval1,err_k_eff
       enddo
       !------------------------------------------
       !Normalized flux --------------------------
             moy = 0.0
             k1 = 1
             j = 1
             do k0 = 1,ng
                 do i=k1,(k1-1)+totNFM
                    moy(k0) = moy(k0) + phi_li1(i)
                 end do
                 do i=k1,(k1-1)+totNFM
                    phi(j) = phi_li1(i)/(moy(k0)/totNFM)
                    j = j + 1
                 end do
                 k1 = k1 + totNFM*(Norder+1)/2
             enddo
       ! Format
       2000 format(t3,"Iteration",i4,":",5x,"===>",5x,"keff =",F9.6,5x,&
                  "===>",5x,"res =",e10.3)
    end subroutine aleig
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
subroutine Matrix_A(Delta,SigR,SigS,fmmid,BC,ng,Nmat,Norder,order,totNFM,ddim,A)
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,totNFM,ddim,Norder,order
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=8), dimension(Nmat,Norder+1,ng), intent(in) :: SigR
       real(kind=8), dimension(totNFM), intent(in) :: Delta
       real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
       CHARACTER(50), intent(in) :: BC
       real(kind=8), dimension(ddim,ddim), intent(out) :: A
       real(kind=8), dimension(totNFM*(Norder+1)/2,totNFM*(Norder+1)/2,ng,ng) :: AA 
       real(kind=8) :: betap=0,betam=0,denom1,denom2,a1,a2,a3,a4,m,pnmars
       real(kind=8), dimension(totNFM,(Norder+1)/2,ng) :: b
       real(kind=8), dimension(totNFM-1,(Norder+1)/2,ng) :: am,ap
       real(kind=8), dimension(Norder) :: be
       real(kind=8), dimension(Norder) :: al
       integer(kind=4) :: i,j,k,k0,k1,k2,k3,k4,k5,n,som,l
       open (100,file='app/Output/matrixa.h')
       b=0.0;am=0.0;ap=0.0;A=0.0;AA=0.0
       !------------------------------------------------------
       ! Prepare coeff beta et alpha
       do n = 1,Norder
          be(n) = float(n)/float(2*n+1)
          al(n) = float(n+1)/float(2*n+1)
       enddo
       print*, be(1)
       !------------------------------------------------------
       ! Select a boundary conditions
       if  ( BC == 'reflective' ) then
           betam = 1.; betap = 1.
       elseif  ( BC == 'vacuum' ) then
           betam = 0.; betap = 0.
       elseif  ( BC == 'vacuum_reflective' ) then
           betam = 0.; betap = 1.
       elseif  ( BC == 'reflective_vacuum' ) then
           betam = 1.; betap = 0.
       endif
       !------------------------------------------------------

       do  l = 1, ng
       n = 2
       do j = 1,(Norder+1)/2
       ! Prepare ak,k-1 ( k= 2 ----> I )
          do i  = 1,totNFM-1
             a1 = Delta(i)*SigR(fmmid(i),n,l)
             a2 = Delta(i+1)*SigR(fmmid(i+1),n,l)
             am(i,j,l)  = -2.0/(a1+a2)  
           enddo
       ! Prepare ak,k+1 ( k= 1 ----> I-1 )
          do i  = 2,totNFM
             a3 = Delta(i-1)*SigR(fmmid(i-1),n,l)
             a4 = Delta(i)*SigR(fmmid(i),n,l)
             ap(i-1,j,l)  = -2.0/(a3+a4)                                         
          enddo
       n = n + 2 
       enddo
       ! Prepare bk,k  ( k= 1 ----> I)
       n = 1
       do  j = 1, (Norder+1)/2
           do i = 1, totNFM
              b(i,j,l)   =  SigR(fmmid(i),n,l)*Delta(i) 
           enddo
       n = n + 2
       enddo
       enddo
      !------------------------------------------------------
       do l = 1,ng
       k = 2; k1=0
       do j = 1, (Norder+1)/2
          n = 1
          do i = k,totNFM*j
             if (j == 1) then
                AA(i-1,i,l,l) = be(j)*am(n,j,l)
                AA(i,i-1,l,l) = be(j)*ap(n,j,l)
             else
                AA(i-1,i,l,l) =  al(k1-1)*be(k1)*am(n,j-1,l)+&
                                  al(k1)*be(k1+1)*am(n,j,l)
                AA(i,i-1,l,l) =  al(k1-1)*be(k1)*ap(n,j-1,l)+&
                                  al(k1)*be(k1+1)*ap(n,j,l)
             endif
             n = n + 1
          enddo
       k1 = k1 + 2
       k = k + totNFM 
       enddo
       ! Prepare Elements in Diagonaux ( ak,k; k= 2 ---> I-1 )
       !------------------------------------------------------
       k = 2;k1=0
       do j = 1,(Norder+1)/2
             n = 1
          do i=k,totNFM*j-1
             if (j == 1) then
             AA(i,i,l,l)  =  b(n+1,j,l)-be(j)*(am(n,j,l)+ap(n,j,l))
             else 
             AA(i,i,l,l)  =  b(n+1,j,l)-al(k1-1)*be(k1)*(am(n,j-1,l)+ap(n,j-1,l))&
                                       -al(k1)*be(k1+1)*(am(n,j,l)+ap(n,j,l))
             endif
             n = n + 1
          enddo  
       k = k + totNFM
       k1 = k1 + 2
       enddo 
       !------------------------------------------------------
       if (Norder > 1) then
          k = 2; k1 = 1; m = al(1); k2 = 2
          do j = 1,(Norder+1)/2-1
             n = 1
             do i=k,totNFM*j-1
                AA(totNFM+i,i,l,l)  =  -be(k1)*be(k1+1)*(am(n,j,l)+ap(n,j,l))
                AA(i,totNFM+i,l,l)  =  -m*(am(n,j,l)+ap(n,j,l))
             n = n + 1
             enddo  
          k= k + totNFM
          m = al(k2)*al(k2+1)
          k1 = k1 + 2
          k2 = k2 + 2
          enddo 
       endif
       !------------------------------------------------------
       if (Norder > 1) then
          k = 2; k1 = 1; m = al(1); k2 = 2
          do j = 1, (Norder+1)/2-1
             n = 1
             do i = k,totNFM*j  
                ! lower 
                AA(totNFM+i-1,i,l,l)  = be(k1)*be(k1+1)*am(n,j,l)
                AA(totNFM+i,i-1,l,l)  = be(k1)*be(k1+1)*ap(n,j,l) 
                ! upper
                AA(i-1,i+totNFM,l,l)  = m*am(n,j,l) 
                AA(i,i-1 +totNFM,l,l) = m*ap(n,j,l)
                n = n + 1
             enddo
          k = k + totNFM 
          m = al(k2)*al(k2+1)
          k1 = k1 + 2
          k2 = k2 + 2
          enddo
       endif
       !------------------------------------------------------
       ! Boundary condition
       
       denom1 = ((1.d0 +betam)/(1.d0 -betam))
       denom2 = ((1.d0 +betap)/(1.d0 -betap))
       k1 = 2; k = totNFM 
       do j = 1, (Norder+1)/2
          if (j == 1) then
              if (Norder == 1) then
              AA(1,1,l,l) = b(1,j,l)+denom1*sqrt(3.0)/3.-be(1)*am(1,j,l)
              AA(k,k,l,l) = b(totNFM,1,l)+denom2*sqrt(3.0)/3.-be(1)*ap(totNFM-1,j,l)
              else
              AA(1,1,l,l) = b(1,j,l)+denom1*pnmars(1,0)-be(1)*am(1,j,l)
              AA(k,k,l,l) = b(totNFM,1,l)+denom2*pnmars(1,0)-be(1)*ap(totNFM-1,j,l)
              endif
          else 
              AA(k+1,k+1,l,l) = b(1,j,l)+denom1*(be(k1)*pnmars(k1-1,k1)+&
                           al(k1)*pnmars(k1+1,k1))-al(k1-1)*be(k1)*am(1,j-1,l)&
                           -al(k1)*be(k1+1)*am(1,j,l)
              AA(k+totNFM,k+totNFM,l,l) = b(totNFM,j,l)+denom2*(be(k1)*pnmars(k1-1,k1)+&
                                   al(k1)*pnmars(k1+1,k1))-al(k1-1)*be(k1)*&
                                ap(totNFM-1,j-1,l)-al(k1)*be(k1+1)*ap(totNFM-1,j,l)
          k1 = k1 + 2
          k = k + totNFM
          endif
       enddo
       !------------------------------------------------------
       if (Norder > 1) then
       k1 = 4; k2=3; k=totNFM ; k3=5
          do j = 1, (Norder+1)/2-1
                if (j==1) then
                AA(k+1,k-totNFM+1,l,l) = denom1*(be(2)*pnmars(1,0)+al(2)*pnmars(3,0))&
                                         -be(1)*be(2)*am(1,j,l)
                AA(k+totNFM,k,l,l) = denom2*(be(2)*pnmars(1,0)+al(2)*pnmars(3,0))&
                                     -be(1)*be(2)*ap(totNFM-1,j,l)
                else
                AA(k+1,k-totNFM+1,l,l) = denom1*(al(k1)*pnmars(k2,k1-2)+be(k1)*&
                                         pnmars(k3,k1-2))-be(k1-1)*be(k1)*am(1,j,l)
                AA(k+totNFM,k,l,l) = denom2*(al(k1)*pnmars(k2,k1-2)+be(k1)*&
                                     pnmars(k3,k1-2))-be(k1-1)*be(k1)*ap(totNFM-1,j,l)
                   if (mod(j,2) == 0) then 
                       k2 = k2 + 4
                   else
                       k3 = k3 + 2
                   endif
                k1 = k1 + 2
                endif
          k = k + totNFM 
          enddo
       endif
       !------------------------------------------------------
       if (Norder > 1) then 
       k1 = 2;k=totNFM
          do j = 1, (Norder+1)/2-1
                if (j==1) then
                AA(k-totNFM+1,k+1,l,l) = denom1*pnmars(1,2)-al(1)*am(1,1,l)                
                AA(k,k+totNFM,l,l) =  denom2*pnmars(1,2)-al(1)*ap(totNFM-1,1,l)
                else
                AA(k-totNFM+1,k+1,l,l) = denom1*(al(k1)*pnmars(k1+1,k1+2)+&
                                         be(k1)*pnmars(k1-1,k1+2))-&
                                         al(k1)*al(k1+1)*am(1,j,l)               
                AA(k,k+totNFM,l,l) =  denom2*(al(k1)*pnmars(k1+1,k1+2)+&
                                      be(k1)*pnmars(k1-1,k1+2))-&
                                      al(k1)*al(k1+1)*ap(totNFM-1,j,l)
                k1 = k1 + 2
                endif
          k = k + totNFM
          enddo
       endif
       !------------------------------------------------------
       if (Norder > 3) then
          !**** Part 1 ****
          k=2;k1=4;k2=(Norder+1)/2-2;k3=6;k5=3
          do j = 1,(Norder+1)/2-2
             som = 0;k0=2;k4=totNFM
             do i=1,k2
                if (i==1) then
                    AA(1,k*totNFM+1,l,l) = denom1*pnmars(1,k1)
                    AA(totNFM,(k+1)*totNFM,l,l) = denom2*pnmars(1,k1)
                    k1 = k1 + 2
                    k = k + 1
                else
                    AA(k4+1,k5*totNFM+1,l,l) = denom1*(al(k0)*pnmars(k0+1,k3)+&
                                          be(k0)*pnmars(k0-1,k3))
                    AA(totNFM+k4,(k5+1)*totNFM,l,l) = denom2*(al(k0)*pnmars(k0+1,k3)+&
                                                 be(k0)*pnmars(k0-1,k3))
                    k0 = k0 + 2
                    som = som + 2
                    k3 = k3 + 2
                    k5 = k5 + 1
                    k4 = k4 + totNFM
                endif
             enddo
             k2 = k2 - 1
             k3 = k3 - som + 2
             k5 = k5 - i + 3
          enddo 
          !**** Part 2 ****
          k0=4;k4=2;k2=(Norder+1)/2-2
          do j=1,(Norder+1)/2-2
             som = 0;k1=0;k5=1
             do i=1,k2
                    AA(k4*totNFM+1,k5,l,l) = denom1*(al(k0)*pnmars(k0+1,k1)+&
                                      be(k0)*pnmars(k0-1,k1))
                    AA((k4+1)*totNFM,(k5-1)+totNFM,l,l) = denom2*(al(k0)*pnmars(k0+1,k1)+&
                                       be(k0)*pnmars(k0-1,k1))
                    som = som + 2
                    k0 = k0 + 2
                    k1 = k1 + 2
                    k4 = k4 + 1
                    k5 = k5 + totNFM
             enddo
             k2 = k2 - 1
             k0 = k0 - som + 2
             k4 = k4 - i + 2
          enddo
       endif
       enddo 
       !------------------------------------------------------

       n=1
       do i = 1,ng
       A(n:n-1+totNFM*(Norder+1)/2,n:n-1+totNFM*(Norder+1)/2) = AA(:,:,i,i) 
       n = n + totNFM*(Norder+1)/2
       enddo
       do i = 1,ddim
       write(100,'(10000f12.8)') (A(j,i), j=1,ddim)   
       enddo
end subroutine Matrix_A
!BL
       function pnmars(l,m)
       implicit none
       integer(kind=4), intent(in) :: l,m
!      return the Marshak boundary coefficients in slab geometry. These
!      coefficients are specific to the left boundary.
!      (C) 2008 Alain Hebert, Ecole Polytechnique de Montreal
       real(kind=8), dimension(16) :: zgksi,wgksi
       real(kind=8) :: pnmars,plgndr,som
       integer(kind=4) :: i
       data zgksi /0.00529953837, 0.0277124941, 0.0671843886, 0.122297794,&
                   0.191061884, 0.270991623, 0.359198213, 0.452493757,&
                   0.547506273, 0.640801787, 0.729008377, 0.808938146,&
                   0.877702236,  0.932815611,  0.972287536,  0.994700432/
       data wgksi /0.01357623,  0.0311267618,  0.047579255,  0.0623144843,&
                   0.0747979954, 0.0845782608, 0.0913017094,  0.0947253034,&
                   0.0947253034,  0.0913017094, 0.0845782608, 0.0747979954,&
                   0.0623144843,  0.047579255,  0.0311267618, 0.01357623/
       if ( mod(l,2) == 0 ) then
       write (*,'(a)') 'error(odd first index expected)'
       stop
       endif
       som=0.0
       do i = 1,16
          som = som + wgksi(i)*plgndr(l,0,zgksi(i))*plgndr(m,0,zgksi(i))
       enddo
       pnmars=(2*m+1)*som 
       return
       end function pnmars
!BL
       function plgndr(l,m,x)
       implicit none
       integer(kind=4), intent(in) :: l,m
       real(kind=8), intent(in) :: x
       real(kind=8) :: plgndr,pmm,pmmp1,pll=0.
       integer(kind=4) :: ll
!      return the Ferrer definition of the associated Legendre function.
!      (C) 2008 Alain Hebert, Ecole Polytechnique de Montreal
       if (m < 0) then
          write (*,'(a)')  'error( bad arguments in plgndr 1)'
          stop
       elseif (m > l) then
          write (*,'(a)')  'error( bad arguments in plgndr 2)'
          stop
       elseif (abs(x) > 1.) then
          write (*,'(a)')  'error( bad arguments in plgndr 3)'
          stop
       endif
       pmm=1. 
       if (l == m) then
          plgndr=pmm 
       else
          pmmp1=(2*m+1)*x*pmm 
          if (l == m+1) then
             plgndr=pmmp1 
          else
          do ll=m+2,l
             pll = ((2*ll-1)*x*pmmp1-(ll+m-1)*pmm)/(ll-m) 
             pmm=pmmp1 ; pmmp1=pll 
          enddo
          plgndr=pll 
          endif
       endif
       end function  plgndr
!BL9
subroutine Matrix_F(NusigF,Chi,F,fmmid,Delta,ng,Nmat,totNFM,Norder,ddim)
       implicit none
       integer(kind=4), intent(in) :: ng,Nmat,totNFM,Norder,ddim
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid
       real(kind=8), dimension(Nmat,ng), intent(in) :: NusigF,Chi
       real(kind=8), dimension(totNFM), intent(in) :: Delta
       real(kind=8), dimension(ddim,ddim), intent(out) :: F
       integer(kind=4) :: i,j,k0,k1,k2,k3,k
       open (101,file='app/Output/matrixf.h')
            F(:,:) = 0.0
            k0 = 1
            do   k   = 1, ng
                 j   = 1
                 do  i  = k0,(k0-1)+totNFM
                     F(i,i) = Delta(j)*Chi(fmmid(j),k)*NusigF(fmmid(j),k)
                     j = j + 1
                 end do
            k0 = totNFM*(Norder+1)/2 + k0
            enddo

            k2  = 0
            do  k1  = 1,ng
                k0  = 1
                k3  = 0
                    do while (k0<k1) 
                       do  i  = 1,totNFM
                       F(i+k3,i+(k1-1)*totNFM*(Norder+1)/2) = Delta(i)*Chi(fmmid(i),k0)*NusigF(fmmid(i),k1)
                       F(i+(k1-1)*totNFM*(Norder+1)/2,i+k3) = Delta(i)*Chi(fmmid(i),k1)*NusigF(fmmid(i),k0)
                       enddo
                       k0 = k0 + 1
                       k3 = k3 + totNFM
                    enddo
            enddo
       do i = 1,ddim 
       write(101,'(10000f12.8)') (F(i,j), j=1,ddim) 
       enddo
end subroutine Matrix_F
!BL
    subroutine Matrix_L(SigS,fmmid,L,ng,totNFM,Nmat,order,Norder,ddim)
       implicit none
       integer(kind=4), intent(in) :: ng,totNFM,Nmat,order,Norder,ddim
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid 
       real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
       real(kind=8), dimension(ddim,ddim), intent(out) :: L
       integer(kind=4) :: i,j,k0,k1,k2,k3
       open (102,file='app/Output/matrixl.h')
                L(:,:) = 0.0
                k2  = 0
            do  k1  = 1,ng
                k0  = 1
                k3  = 0
                    do while (k0<k1) 
                             do  i  = 1,totNFM
                             L(i+(k1-1)*totNFM*(Norder+1)/2,i+k3) = SigS(fmmid(i),1,k0,k1) ! 1 --> 2 SigS(1,1,2)
                             enddo
                             k0 = k0 + 1
                             k3 = k3 + totNFM
                    enddo
            enddo
       do i = 1,ddim 
       write(102,'(10000f12.8)') (L(i,j), j=1,ddim) 
       enddo
    end subroutine 
!BL    
    subroutine Matrix_U(SigS,fmmid,U,ng,totNFM,Nmat,order,Norder,ddim)
       implicit none
       integer(kind=4), intent(in) :: ng,totNFM,Nmat,order,Norder,ddim
       integer(kind=4), dimension(totNFM), intent(in) :: fmmid 
       real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
       real(kind=8), dimension(ddim,ddim), intent(out) :: U
       integer(kind=4) :: i,j,k0,k1,k2,k3
       open (103,file='app/Output/matrixu.h')
                U(:,:) = 0.0
                k2  = 0
            do  k1  = 1,ng
                k0  = 1
                k3  = 0
                    do while (k0<k1) 
                             do  i  = 1,totNFM
                             U(i+k3,i+(k1-1)*totNFM*(Norder+1)/2) = SigS(fmmid(i),1,k1,k0) ! 1 --> 2 SigS(1,1,2)
                             enddo
                             k0 = k0 + 1
                             k3 = k3 + totNFM
                    enddo
            enddo 
       do i = 1,ddim 
       write(103,'(10000f12.8)') (U(i,j), j=1,ddim) 
       enddo
    end subroutine Matrix_U
!BL
subroutine RemovalXS(ng,Nmat,order,Norder,SigT,SigS,SigR)
       implicit none
       integer(kind=4), intent(in) ::  ng,Nmat,order,Norder
       real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
       real(kind=8), dimension(Nmat,ng), intent(in) :: SigT
       real(kind=8), dimension(Nmat,Norder+1,ng), intent(out) :: SigR
       integer(kind=4) :: i,k0,k1,k4,k
       do k1 = 1,order
            do   k   = 1, ng
                 do  i  = 1,Nmat
                     SigR(i,k1,k) = SigT(i,k) - SigS(i,k1,k,k)
                 end do
            enddo
       enddo

       do k1 = order+1,Norder+1
            do   k   = 1, ng
                 do  i  = 1,Nmat
                     SigR(i,k1,k) = SigT(i,k) 
                 end do
            enddo
       enddo
       !write(*,'(2f8.4)') transpose(SigR(1,1,:,:))
end subroutine RemovalXS
!BL10
subroutine flux_guess(totNFM,ng,ddim,Nmat,Norder,nregion,NusigF,dcell,phi_li)
       ! Guess scalar flux
       implicit none
       integer(kind=4), intent(in) ::  totNFM,ng,ddim,Nmat,nregion,Norder
       real(kind=8), dimension(Nmat,ng), intent(in) :: NusigF
       real(kind=8), dimension(nregion), intent(in) :: dcell
       real(kind=8), dimension(ddim), intent(out) :: phi_li
       integer(kind=4) :: i,k,k0,k1,k2
       real(kind=8), dimension(ng*nregion) :: a10,a11
       real(kind=8) :: const
       i = 1
       phi_li  = 0.0
       do k1=1,ng
          do K2=1,nregion
             a10(i) = NusigF(k2,k1)
             a11(i) = dcell(k2)
             i = i + 1
          enddo
       enddo 
       const = 1.0/dot_product(a10,a11)

       k0 = 1
       do   k   = 1, ng
            do  i  = k0,(k0-1)+totNFM
                phi_li(i)  = const
            end do
            k0 = totNFM*(Norder+1)/2 + k0
       enddo
       !write(*,'(10f8.4)') phi_li
end subroutine flux_guess
!BL13
subroutine Output(start,BC,tm,k_eff,SigT,NusigF,SigS,Chi,dcell,phi,eps,totNFM,dim,&
                  ng,Nmat,order,nregion,norder,it1,it2)
        implicit none
        integer(kind=4), intent(in) :: ng,dim,totNFM,Nmat,order,nregion,norder,it1,it2
        real(kind=8), dimension(Nmat,ng), intent(in) :: SigT,NusigF,Chi
        real(kind=8), dimension(Nmat,order,ng,ng), intent(in) :: SigS
        real(kind=8), dimension(dim), intent(in) :: phi
        real(kind=8), dimension(nregion), intent(in) :: dcell
        CHARACTER(50), intent(in) :: start,BC,tm
        real(kind=8), intent(in) :: eps,k_eff
        real(kind=8) :: pnmars !pnmars(l,m)
        integer(kind=4) :: i,j
        open (100,file='app/Output/OUTPUT_PN.TXT')
        write (100, FMT=* ) '********************************************************************************'
        write (100, FMT=* ) 'ERSN, UNIVERSITY ABDELMALEK ESSAADI FACULTY OF SCIENCES - TETOUAN, MOROCCO'
        write (100, FMT=* ) 'CODE  DEVELOPED  BY  MOHAMED  LAHDOUR,  PHD  STUDENT'
        write (100, FMT=* ) 'NTP-ERSN:        PN  SPHERICAL HARMONIC METHOD'
        write (100, FMT=* ) 'VERSION NUMBER:  1.2'
        write (100, FMT=* ) 'VERSION DATE:    10  OTOBER  2018'
        write (100,3010) 'RAN ON:          ', start,'(H:M:S)'
        write (100, FMT=* ) '********************************************************************************'
        write (100, FMT=* ) '           ----------------------------------------------------------' 
        write (100, FMT=* ) '                     INPUT  PARAMETER - VALUES  FROM  INPUT'              
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'ENERGY GROUPS NUMBER:                  ',ng
        write (100, FMT=* ) 'REGIONS NUMBER:                        ',nregion
        write (100, FMT=* ) 'MATERIALS NUMBER:                      ',Nmat
        write (100,3040)    'SIZE OF EACH REGION [CM]:              ',dcell       
        write (100, FMT=* ) 'L-ORDER LEGENDRE POLYNOMIAL:           ',order-1 
        write (100, FMT=* ) 'N-ORDER LEGENDRE POLYNOMIAL:           ',norder
        write (100, FMT=* ) 'TOTAL NUMBER OF FINE MESHES:           ',totNFM
        write (100,3050)    'CONVERGENCE CRITERION of KEFF AND FLUX:',eps
        write (100, FMT=* ) ''
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) '                      CALCULATION  RUN-TIME  PARANETERS  PN' 
        write (100, FMT=* ) '           ----------------------------------------------------------'
        write (100, FMT=* ) ''
        write (100, FMT=* ) 'MARSHAK COEFFICIENTS PNMARS(l,m) IN SLAB GEOMETRY: '
        write (100, FMT=* ) ''
        write (100,2020) '       l ',('     m =',i,i=0,norder,2)
        write (100, FMT=* ) ''
        print*,norder
        do i=1,norder,2
        write(100,2010) i,(pnmars(i,j), j=0,norder,2)
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
        2010 format(1x,1p,i11,5x,200e16.5)
        2020 format(1x,A12,1x,300(A14,i2))
        3000 format(1x,A14,2x,300(A14,i2))  
        3010 format(1x,A17,A22,A10)
        3020 format(1x,A26,4x,i10)
        3040 format(1x,A39,2x,200F10.5)
        3050 format(1x,1p,A39,4x,e8.1)
        3060 format(1x,1p,i11,5x,e16.5,e16.5)
        3070 format(1x,A18,i4)
        3080 format(1x,1p,i11,5x,e16.5,e16.5,e16.5,e16.5,e16.5)
        3090 format(1x,A26,6x,f8.6)
        4000 format(1x,A26,4x,A10,A10)
        close(100)
end subroutine Output
!BL14
subroutine plot_flux(Delta,flux,totNFM,dim)
        implicit none
        integer(kind=4), intent(in) :: dim,totNFM
        real(kind=8), dimension(dim), intent(in) :: flux
        real(kind=8), dimension(totNFM), intent(in) :: Delta
        real(kind=8) :: som
        integer(kind=4) :: i,j
        open (10,file='app/Output/flux_pn.h')
        som = Delta(1)
        do i=1,totNFM
        write(10,'(10000f12.8)') som,(flux(i+j), j=0,dim-1,totNFM)  
        som = som + Delta(i)
        enddo
        close(10)
end subroutine plot_flux
!BL15
subroutine title1(CM)   
       CHARACTER(80), intent(in) :: CM    
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
       write(*,FMT=*)'                                                   Version Number: 1.2 '
       write(*,FMT=*)'     Copyright:      2015-2018 FS-Tetouan University Abdelmalk Essaadi '
       write ( *, FMT=* ) '. '
       write ( *, FMT=* ) '   FORTRAN90 version'  
       write ( *, FMT=* ) '   The Spherical Harmonic Method Pn'  
       write ( *, FMT=* ) '   Calculation of 1D Discrete Moment and Space Domain '
       write ( *, FMT=* ) '   Slab 1D geometry'
       write ( *, FMT=* ) '   Solution of Linear Algebraic Equations by :'
       write ( *, FMT=* ) '   ',CM
       write ( *, FMT=* ) '. '
end subroutine title1
!BL16
subroutine title2()
       write ( *, FMT=* )' ************************************************************************'
       write ( *, FMT=* )'                               Finished'                             
       write ( *, FMT=* )' ************************************************************************'  
end subroutine title2
!BL27
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
!BL 
      subroutine LU_decomp(n,A,L,U)
      implicit none
      integer(kind=4), intent(in) :: n
      real(kind=8), dimension(n,n), intent(in) :: A
      real(kind=8), dimension(n,n), intent(out) :: L,U
      integer(kind=4) :: i,j,k
      do k=1,n
         L(k,k) = A(k,k) - sum(L(k,1:k-1)*U(1:k-1,k))
         do j = k,n
            U(k,j) = (A(k,j)-sum(L(k,1:k-1)*U(1:k-1,j)))/L(k,k)
         enddo
         do i = k+1,n
            L(i,k) = (A(i,k)-sum(L(i,1:k-1)*U(1:k-1,k)))/U(k,k)
         enddo
       enddo
      end subroutine LU_decomp
!BL
      subroutine ForSub(n,L,B,Y)
      implicit none
      integer(kind=4), intent(in) :: n
      real(kind=8), dimension(n), intent(in) :: B
      real(kind=8), dimension(n,n), intent(in) :: L
      real(kind=8), dimension(n), intent(out) :: Y
      integer(kind=4)  :: i
      Y(1) = B(1)/L(1,1)
      do i = 2,n
          Y(i) = (B(i)-sum(L(i,1:i-1)*Y(1:i-1)))/L(i,i)
      enddo 
      end subroutine ForSub
!BL
      subroutine BackSub(n,U,Y,X)
      implicit none
      integer(kind=4), intent(in) :: n
      real(kind=8), dimension(n), intent(in) :: Y
      real(kind=8), dimension(n,n), intent(in) :: U
      real(kind=8), dimension(n), intent(out) :: X
      integer(kind=4)  :: i
      X(n) = Y(n)/U(n,n)
      do i = n-1,1,-1
          X(i) = (Y(i)-sum(U(i,i+1:n)*X(i+1:n)))/U(i,i)
      enddo 
      end subroutine BackSub































