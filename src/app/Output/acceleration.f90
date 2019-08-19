 PROGRAM  ITVMET
      PARAMETER (N=3)
      INTEGER::I,J
      REAL::A(10,10),A1(10,10),A2(10,10),B(10),B1(10),B2(10)
      REAL::X0(10),X01(10),X02(10),TOL,W
      OPEN(1,’INPUT.DAT’)
      OPEN(2,’OUTPUT.DAT’)
      READ(1,*)((A(I,J),J=1,N),I=1,N)
      READ(1,*)(B(I),I=1,N)
      READ(1,*)(X0(I),I=1,N)
      READ(1,*)TOL,W
      DO I=1,N
      B1(I)=B(I)
      B2(I)=B(I)
      X01(I)=X0(I)
      X02(I)=X0(I)
      DO J=1,N
      A1(I,J)=A(I,J)
      A2(I,J)=A(I,J)
      END DO
      END DO
      CALL JM(A,B,X0,N,TOL)
      CALL GSM(A1,B1,X01,N,TOL)
      CALL SOR(A2,B2,X02,N,TOL,W)
      END PROGRAM

      SUBROUTINE JM(A,B,X0,N,TOL)
      REAL::A(10,10),B(10),X0(10),X(10),NORM,SUM1
      INTEGER::K=1
      WRITE (2,*)’RESULT FOR JACOBI METHOD’
  10  DO I=1,N
      SUM1=0.0
      DO J=1,N
      IF (J.NE.I) SUM1=SUM1+A(I,J)*X0(J)
      END DO
      X(I)=(B(I)-SUM1)/A(I,I)
      END DO
      WRITE (2,12)K,(X(I),I=1,N)
  12  FORMAT(2X,I3,3(2X,F9.6))
      K=K+1
      NORM=ABS(X(1)-X0(1))
      DO I=2,N
      IF (ABS(X(I)-X0(I)).GT.NORM) NORM=ABS(X(I)-X0(I))
      END DO
      IF (NORM.LT.TOL) GOTO 11
      DO I=1,N
      X0(I)=X(I)
      END DO
      GO TO 10
   11 END SUBROUTINE

      SUBROUTINE GSM(A1,B1,X01,N,TOL)
      REAL::A1(10,10),B1(10),X01(10),X(10),NORM,SUM1,SUM2
      INTEGER::K=1
      WRITE (2,*)’RESULT FOR GAUSS-SEIDEL METHOD’
  11  DO I=1,N
      SUM1=0.0
      SUM2=0.0
      DO J=1,N
      IF (J.LT.I) SUM1=SUM1+A1(I,J)*X(J)
      IF (J.GT.I) SUM2=SUM2+A1(I,J)*X01(J)
      END DO
      X(I)=(B1(I)-SUM1-SUM2)/A1(I,I)
      END DO
      WRITE (2,20)K,(X(I),I=1,N)
  20  FORMAT(2X,I3,3(2X,F9.6))
      K=K+1
      NORM=ABS(X(1)-X01(1))
      DO I=2,N
      IF (ABS(X(I)-X01(I)).GT.NORM) NORM=ABS(X(I)-X01(I))
      END DO
      IF (NORM.LT.TOL) GO TO 12
      DO I=1,N
      X01(I)=X(I)
      END DO
      GO TO 11
   12 END SUBROUTINE

      SUBROUTINE SOR(A2,B2,X02,N,TOL,W)
      REAL::A2(10,10),B2(10),X02(10),X(10),NORM,SUM1,SUM2,W
      INTEGER::K=1
      WRITE(2,*)’RESULT FOR SOR METHOD’
  13  DO I=1,N
      SUM1=0.0
      SUM2=0.0
      DO J=1,N
      IF (J.LT.I) SUM1=SUM1+A2(I,J)*X(J)
      IF (J.GT.I) SUM2=SUM2+A2(I,J)*X02(J)
      END DO
      X(I)=(1.0-W)*X02(I)+(W*(B2(I)-SUM1-SUM2))/A2(I,I)
      END DO
      WRITE(2,30)K,(X(I),I=1,N)
  30  FORMAT(2X,I3,3(2X,F9.6))
      K=K+1
      NORM=ABS(X(1)-X02(1))
      DO I=2,N
      IF(ABS(X(I)-X02(I)).GT.NORM) NORM=ABS(X(I)-X02(I))
      END DO
      IF (NORM.LT.TOL) GO TO 14
      DO I=1,N
      X02(I)=X(I)
      END DO
      GO TO 13
  14  END SUBROUTINE
