      SUBROUTINE SNQU02(NLF,xi,mu,eta,we)
! 
!-----------------------------------------------------------------------
! 
!Purpose:
! set the level-symmetric (type 2) quadratures.
! 
!Copyright:
! Copyright (C) 2005 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version
! 
!Author(s): A. Hebert
! 
!Parameters: input
! NLF     order of the SN approximation (even number).
! 
!Parameters: output
! JOP     number of base points per axial level in one octant.
! U       base points in $\xi$ of the axial quadrature. Used with
!         zero-weight points.
! W       weights for the axial quadrature in $\xi$.
! TPQ     base points in $\xi$ of the 2D SN quadrature.
! UPQ     base points in $\mu$ of the 2D SN quadrature.
! VPQ     base points in $\eta$ of the 2D SN quadrature.
! WPQ     weights of the 2D SN quadrature.
! 
!-----------------------------------------------------------------------
! 
      !IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!----
!  SUBROUTINE ARGUMENTS
!----
      impl
icit none
      integer(kind=8),  intent(in) :: NLF
      real(kind=8), dimension(NLF*(NLF/2+1)/4), intent(out) :: xi,mu
      real(kind=8), dimension(NLF*(NLF/2+1)/4), intent(out) :: eta
      real(kind=8), dimension(NLF*(NLF/2+1)/4), intent(out) :: we
      real(kind=8), dimension(NLF/2) :: U,W
      integer(kind=8), dimension(NLF/2) :: JOP
!----
!  LOCAL VARIABLES
!----
      REAL(kind=8), PARAMETER :: PI = 3.141592654
      integer(kind=8) ::  M2,NPQ,I,IPR,INMAX,IP,IQ,IS,NW0,KK,LL,II
      integer(kind=8) ::  IPK,IPL,IPQ,IW,NEQ,NW
      real(kind=8) :: ZMU1,ZMU2
      INTEGER(kind=8), PARAMETER :: MAXNLF=24,MAXEQ=64
      INTEGER(kind=8), PARAMETER :: MAXW = int(1+(MAXNLF*&
                                  (MAXNLF+8)-1)/float(48))
      INTEGER(kind=8), dimension(MAXNLF*(MAXNLF/2+1)/4) :: INWEI
      DOUBLE PRECISION :: WSUM2, WSUM, REF,ZETA,ZMU,ZMAT(MAXEQ,MAXW+1),UD(MAXW),WEI(MAXW)

      
!----
!  SET THE UNIQUE QUADRATURE VALUES.
!----
      IF(NLF.GT.MAXNLF) print*,'SNQU02: MAXNLF OVERFLOW.'
      M2=NLF/2
      NPQ=M2*(M2+1)/2 
      ZMU1=1.0D0/(3.0D0*DBLE(NLF-1)) ! The first base point 
      NW=1+(NLF*(NLF+8)-1)/48 ! Number of distinct weights per octant 
      IF(NLF.EQ.2) THEN
         ZMU1=0.33333333
      ELSE IF(NLF.EQ.4) THEN
         ZMU1=0.12251480
      ELSE IF(NLF.EQ.6) THEN
         ZMU1=0.07109447
      ELSE IF(NLF.EQ.8) THEN
         ZMU1=0.04761903
      ELSE IF(NLF.EQ.10) THEN
         ZMU1=0.03584310
      ELSE IF(NLF.EQ.12) THEN
         ZMU1=0.02796615
      ELSE IF(NLF.EQ.14) THEN
         ZMU1=0.02310250
      ELSE IF(NLF.EQ.16) THEN
         ZMU1=0.01931398
      ELSE IF(NLF.EQ.18) THEN
         ZMU1=0.01692067
      ELSE IF(NLF.EQ.20) THEN
         ZMU1=0.01455253
      ELSE
         print*,'SNQU02: ORDER NOT AVAILABLE FOR N.'
         stop
      ENDIF
      U(1)=REAL(SQRT(ZMU1))
      DO I=2,M2
         ZMU2=ZMU1+2.0D0*DBLE(I-1)*(1.0D0-3.0D0*ZMU1)/DBLE(NLF-2)
         U(I)=REAL(SQRT(ZMU2))
      ENDDO

!----
!  COMPUTE THE POSITION OF WEIGHTS.
!----
      
      IPR=0
      INMAX=0
      DO IP=1,M2
         JOP(IP)=M2-IP+1
         DO IQ=1,JOP(IP)
            IPR=IPR+1
            xi(IPR)=U(IP)
            mu(IPR)=U(M2+2-IP-IQ)
            eta(IPR)=U(IQ)
            IS=MIN(IP,IQ,M2+2-IP-IQ)
            NW0=0
            DO II=1,IS-1
               NW0=NW0+(M2-3*(II-1)+1)/2
            ENDDO
            KK=IP-IS+1
            LL=IQ-IS+1
            IF(KK.EQ.1)THEN
               INWEI(IPR)=NW0+MIN(LL,M2-3*(IS-1)+1-LL)
            ELSEIF(LL.EQ.1)THEN
               INWEI(IPR)=NW0+MIN(KK,M2-3*(IS-1)+1-KK)
            ELSE
               INWEI(IPR)=NW0+MIN(KK,LL)
            ENDIF
            INMAX=MAX(INMAX,INWEI(IPR))
         ENDDO
      ENDDO

      IF(INMAX.NE.NW) print*,'SNQU02: INVALID VALUE OF NW.'
      IF(IPR.NE.NPQ)  print*,'SNQU02: BAD VALUE ON NPQ.'

!----
!  SET THE RECTANGULAR SYSTEM AND SOLVE IT USING THE QR METHOD.
!----
      NEQ=0
      DO IPL=0,NLF,2
        DO IPK=IPL,NLF-IPL,2
          IF(MOD(IPL+IPK,2).EQ.1) CYCLE
          NEQ=NEQ+1
          IF(NEQ.GT.MAXEQ) print*,'SNQU02: MAXEQ OVERFLOW.'
          DO IW=1,NW
             ZMAT(NEQ,IW)=0.0D0
          ENDDO
          DO IPQ=1,NPQ
             ZMU=xi(IPQ)
             ZETA=mu(IPQ)
             IW=INWEI(IPQ)
             ZMAT(NEQ,IW)=ZMAT(NEQ,IW)+(ZMU**IPK)*(ZETA**IPL)
          ENDDO
          REF=1.0D0/DBLE(IPK+IPL+1)
          DO I=1,IPL-1,2
             REF=REF*DBLE(I)/DBLE(IPK+I)
          ENDDO
          ZMAT(NEQ,NW+1)=REF
        ENDDO
      ENDDO
      CALL ALST2F(MAXEQ,NEQ,NW,ZMAT,UD)
      CALL ALST2S(MAXEQ,NEQ,NW,ZMAT,UD,ZMAT(1,NW+1),WEI)


!----
!  SET THE LEVEL-SYMMETRIC QUADRATURES.
!----
      IPQ=0
      WSUM=0.0
      DO IP=1,M2
         WSUM2=0.0D0
         DO IQ=1,JOP(IP)
            IPQ=IPQ+1
            we(IPQ)=REAL(WEI(INWEI(IPQ))*PI/2.0)
            WSUM2=WSUM2+WEI(INWEI(IPQ))
         ENDDO
         W(IP)=REAL(WSUM2)
         WSUM=WSUM+REAL(WSUM2*PI/2.0)
      ENDDO

      RETURN
      END SUBROUTINE

      SUBROUTINE ALST2F(MDIM,M,N,A,UD)
! 
!-----------------------------------------------------------------------
! 
!  TO OBTAIN THE Q*U DECOMPOSITION OF THE MATRIX A USING HOUSEHOLDER
!  TRANSFORMATIONS.
! 
!  INPUT PARAMETERS:
!    MDIM - THE DIMENSIONED COLUMN LENGTH OF A
!    M    - THE NUMBER OF ROWS IN A.
!    N    - THE NUMBER OF COLUMNS IN A. N.LE.M IS ASSUMED.
!    A    - THE MATRIX A.
! 
!  OUTPUT PARAMETER:
!    A  - THE DECOMPOSED MATRIX.
!         LET V(J)(I)=0      FOR I=1,...,J-1
!         AND V(J)(I)=A(I,J) FOR I=J,...,M, THEN
! 
!             Q = PRODUCT(J=1,...,N)(I-BETA(J)*V(J)*V(J)-TRANSPOSE)
! 
!         WHERE BETA(J)=1/(UD(J)*ABS(A(J,J)))
!         A(I,J) FOR I.LT.J GIVES THE OFF-DIAGONAL ELEMENTS OF U.
!    UD - THE DIAGONAL OF U.
! 
!  SCRATCH SPACE ALLOCATED - NONE.
! 
!  ERROR STATES -
! 
!    1 - MDIM.LT.M.
!    2 - N.LT.1.
!    3 - N.GT.M.
!    4 - A IS RANK-DEFICIENT. (RECOVERABLE)
! 
!  P.A. BUSINGER, NUM. MATH. 7, 269-276 (1965).
! 
!-----------------------------------------------------------------------
! 
      integer(kind=8),  intent(in) :: M,N,MDIM
      DOUBLE PRECISION, dimension(N),intent(inout) :: UD
      DOUBLE PRECISION, dimension(MDIM,N),intent(inout) :: A
      DOUBLE PRECISION S,BE,DSQRT
      CHARACTER HSMG*131
      integer(kind=8) ::  I,J,L
! 
! ... CHECK THE INPUT.
! 
      IF(MDIM.LT.M) print*,'ALST2F: MDIM.LT.M'
      IF(N.LT.1) print*,'ALST2F: N.LT.1'
      IF(N.GT.M) THEN
         WRITE(HSMG,'(18HALST2F: N.GT.M (N=,I3,3H M=,I3,2H).)') N,M
         !CALL XABORT(HSMG)
      ENDIF
! 
      DO  L=1,N
      S=0.D0
      DO  I=L,M
      S=S+A(I,L)**2
      ENDDO
      IF(S.EQ. 0.0D0) THEN
         print*,'ALST2F: A IS RANK-DEFICIENT'
      ENDIF
      UD(L)=DSQRT(S)
      IF(A(L,L) .GE. 0.0D0) UD(L)=-UD(L)
      BE=S-A(L,L)*UD(L)
      A(L,L)=A(L,L)-UD(L)
      IF(L.EQ.N) GO TO 50
      DO  J=L+1,N
          S=0.D0
          DO  I=L,M
          S=S+A(I,L)*A(I,J)
          ENDDO
          S=S/BE
          DO  I=L,M
          A(I,J)=A(I,J)-S*A(I,L)
          ENDDO
      ENDDO
      ENDDO
 50   RETURN
      END

      SUBROUTINE ALST2S(MDIM,M,N,A,UD,B,X)
! 
!-----------------------------------------------------------------------
! 
! TO SOLVE THE LEAST SQUARES PROBLEM A*X=B WHEN THE MATRIX A HAS
! ALREADY BEEN DECOMPOSED BY ALST2F.
! 
!  INPUT PARAMETERS:
!    MDIM - THE DIMENSIONED COLUMN LENGTH OF A.
!    M    - THE NUMBER OF ROWS OF A
!    N    - THE NUMBER OF COLUMNS OF A. N.LE.M IS ASSUMED.
!    A    - THE DECOMPOSED MATRIX.
!    UD   - THE DIAGONAL OF U.
!    B    - THE RIGHT-HAND SIDE.
! 
!  OUTPUT PARAMETERS:
!    B - B HAS BEEN CLOBBERED.
!        SQRT(SUM(I=N+1,M)(B(I)**2)) IS THE L2 NORM OF THE RESIDUAL
!        IN THE SOLUTION OF THE EQUATIONS.
!    X - THE SOLUTION VECTORS. X=B IS OK.
! 
!  SCRATCH SPACE ALLOCATED - NONE.
! 
!  ERROR STATES -
! 
!    1 - MDIM.LT.M.
!    2 - N.LT.1.
!    3 - N.GT.M.
!    4 - UD(J)=0 OR A(J,J)=0.
! 
!-----------------------------------------------------------------------
! 
      integer(kind=8),  intent(in) :: M,N,MDIM
      DOUBLE PRECISION, dimension(N),intent(inout) :: UD,X
      DOUBLE PRECISION, dimension(MDIM,N),intent(inout) :: A
      DOUBLE PRECISION, dimension(M),intent(out) :: B
      DOUBLE PRECISION S
      integer(kind=8) ::  I,J,II
! 
! ... CHECK THE INPUT.
! 
      IF(MDIM.LT.M) print*,'ALST2S: MDIM.LT.M'
      IF(N.LT.1) print*,'ALST2S: N.LT.1'
      IF(N.GT.M) print*,'ALST2S: N.GT.M'
! 
! ... APPLY Q-TRANSPOSE TO B.
! 
      DO  J=1,N
      IF((UD(J).EQ. 0.0D0).OR.(A(J,J).EQ. 0.0D0)) THEN
         print*,'ALST2S: UD(J)=0 OR A(J,J)=0'
      ENDIF
      S=0.D0
      DO  I=J,M
      S=S+A(I,J)*B(I)
      ENDDO
      S=S/(UD(J)*A(J,J))
      DO  I=J,M
      B(I)=B(I)+S*A(I,J)
      ENDDO
      ENDDO
! 
! ... BACK-SOLVE THE TRIANGULAR SYSTEM U*X=(Q-TRANSPOSE)*B.
! 
      X(N)=B(N)/UD(N)
      IF(N.EQ.1) GO TO 50
      DO  II=2,N
      I=N+1-II
      S=B(I)
      DO J=I+1,N
         S=S-A(I,J)*X(J)
      ENDDO
      X(I)=S/UD(I)
      ENDDO
! 
 50   RETURN
      END
