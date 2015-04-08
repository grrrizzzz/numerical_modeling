C   IMSL ROUTINE NAME   - FFT2C
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - ELXSI/SINGLE
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - COMPUTE THE FAST FOURIER TRANSFORM OF A
C                           COMPLEX VALUED SEQUENCE OF LENGTH EQUAL TO
C                           A POWER OF TWO
C
C   USAGE               - CALL FFT2C (A,M,IWK)
C
C   ARGUMENTS    A      - COMPLEX VECTOR OF LENGTH N, WHERE N=2**M.
C                           ON INPUT A CONTAINS THE COMPLEX VALUED
C                           SEQUENCE TO BE TRANSFORMED.
C                           ON OUTPUT A IS REPLACED BY THE
C                           FOURIER TRANSFORM.
C                M      - INPUT EXPONENT TO WHICH 2 IS RAISED TO
C                           PRODUCE THE NUMBER OF DATA POINTS, N
C                           (I.E. N = 2**M).
C                IWK    - WORK VECTOR OF LENGTH M+1.
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  FFT2C COMPUTES THE FOURIER TRANSFORM, X, ACCORDING
C                TO THE FOLLOWING FORMULA;
C
C                  X(K+1) = SUM FROM J = 0 TO N-1 OF
C                           A(J+1)*CEXP((0.0,(2.0*PI*J*K)/N))
C                  FOR K=0,1,...,N-1 AND PI=3.1415...
C
C                NOTE THAT X OVERWRITES A ON OUTPUT.
C            2.  FFT2C CAN BE USED TO COMPUTE THE INVERSE FOURIER
C                TRANSFORM, X, ACCORDING TO THE FOLLOWING FORMULA;
C
C                  X(K+1) = (1/N)*SUM FROM J = 0 TO N-1 OF
C                           A(J+1)*CEXP((0.0,(-2.0*PI*J*K)/N))
C                  FOR K=0,1,...,N-1 AND PI=3.1415...
C
C                BY PERFORMING THE FOLLOWING STEPS;
C
C                     DO 10 I=1,N
C                        A(I) = CONJG(A(I))
C                  10 CONTINUE
C                     CALL FFT2C (A,M,IWK)
C                     DO 20 I=1,N
C                        A(I) = CONJG(A(I))/N
C                  20 CONTINUE
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE  FFT2C (A,M,IWK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,IWK(1)
      COMPLEX            A(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ISP,J,JJ,JSP,K,K0,K1,K2,K3,KB,KN,MK,MM,MP,N,
     1                   N4,N8,N2,LM,NN,JK
      REAL               RAD,C1,C2,C3,S1,S2,S3,CK,SK,SQ,A0,A1,A2,A3,
     1                   B0,B1,B2,B3,TWOPI,TEMP,
     2                   ZERO,ONE,Z0(2),Z1(2),Z2(2),Z3(2)
      COMPLEX            ZA0,ZA1,ZA2,ZA3,AK2
      EQUIVALENCE        (ZA0,Z0(1)),(ZA1,Z1(1)),(ZA2,Z2(1)),
     1                   (ZA3,Z3(1)),(A0,Z0(1)),(B0,Z0(2)),(A1,Z1(1)),
     2                   (B1,Z1(2)),(A2,Z2(1)),(B2,Z2(2)),(A3,Z3(1)),
     3                   (B3,Z3(2))
      DATA               SQ/.7071068/,
     1                   SK/.3826834/,
     2                   CK/.9238795/,
     3                   TWOPI/6.283185/
      DATA               ZERO/0.0/,ONE/1.0/
C                                  SQ=SQRT2/2,SK=SIN(PI/8),CK=COS(PI/8)
C                                  TWOPI=2*PI
C                                  FIRST EXECUTABLE STATEMENT
      MP = M+1
      N = 2**M
      IWK(1) = 1
      MM = (M/2)*2
      KN = N+1
C                                  INITIALIZE WORK VECTOR
      DO 5  I=2,MP
         IWK(I) = IWK(I-1)+IWK(I-1)
    5 CONTINUE
      RAD = TWOPI/N
      MK = M - 4
      KB = 1
      IF (MM .EQ. M) GO TO 15
      K2 = KN
      K0 = IWK(MM+1) + KB
   10 K2 = K2 - 1
      K0 = K0 - 1
      AK2 = A(K2)
      A(K2) = A(K0) - AK2
      A(K0) = A(K0) + AK2
      IF (K0 .GT. KB) GO TO 10
   15 C1 = ONE
      S1 = ZERO
      JJ = 0
      K = MM - 1
      J = 4
      IF (K .GE. 1) GO TO 30
      GO TO 70
   20 IF (IWK(J) .GT. JJ) GO TO 25
      JJ = JJ - IWK(J)
      J = J-1
      IF (IWK(J) .GT. JJ) GO TO 25
      JJ = JJ - IWK(J)
      J = J - 1
      K = K + 2
      GO TO 20
   25 JJ = IWK(J) + JJ
      J = 4
   30 ISP = IWK(K)
      IF (JJ .EQ. 0) GO TO 40
C                                  RESET TRIGONOMETRIC PARAMETERS
      C2 = JJ * ISP * RAD
      C1 = COS(C2)
      S1 = SIN(C2)
   35 C2 = C1 * C1 - S1 * S1
      S2 = C1 * (S1 + S1)
      C3 = C2 * C1 - S2 * S1
      S3 = C2 * S1 + S2 * C1
   40 JSP = ISP + KB
C                                  DETERMINE FOURIER COEFFICIENTS
C                                    IN GROUPS OF 4
      DO 50 I=1,ISP
         K0 = JSP - I
         K1 = K0 + ISP
         K2 = K1 + ISP
         K3 = K2 + ISP
         ZA0 = A(K0)
         ZA1 = A(K1)
         ZA2 = A(K2)
         ZA3 = A(K3)
         IF (S1 .EQ. ZERO) GO TO 45
         TEMP = A1
         A1 = A1 * C1 - B1 * S1
         B1 = TEMP * S1 + B1 * C1
         TEMP = A2
         A2 = A2 * C2 - B2 * S2
         B2 = TEMP * S2 + B2 * C2
         TEMP = A3
         A3 = A3 * C3 - B3 * S3
         B3 = TEMP * S3 + B3 * C3
   45    TEMP = A0 + A2
         A2 = A0 - A2
         A0 = TEMP
         TEMP = A1 + A3
         A3 = A1 - A3
         A1 = TEMP
         TEMP = B0 + B2
         B2 = B0 - B2
         B0 = TEMP
         TEMP = B1 + B3
         B3 = B1 - B3
         B1 = TEMP
         A(K0) = CMPLX(A0+A1,B0+B1)
         A(K1) = CMPLX(A0-A1,B0-B1)
         A(K2) = CMPLX(A2-B3,B2+A3)
         A(K3) = CMPLX(A2+B3,B2-A3)
   50 CONTINUE
      IF (K .LE. 1) GO TO 55
      K = K - 2
      GO TO 30
   55 KB = K3 + ISP
C                                  CHECK FOR COMPLETION OF FINAL
C                                    ITERATION
      IF (KN .LE. KB) GO TO 70
      IF (J .NE. 1) GO TO 60
      K = 3
      J = MK
      GO TO 20
   60 J = J - 1
      C2 = C1
      IF (J .NE. 2) GO TO 65
      C1 = C1 * CK + S1 * SK
      S1 = S1 * CK - C2 * SK
      GO TO 35
   65 C1 = (C1 - S1) * SQ
      S1 = (C2 + S1) * SQ
      GO TO 35
   70 CONTINUE
C                                  PERMUTE THE COMPLEX VECTOR IN
C                                    REVERSE BINARY ORDER TO NORMAL
C                                    ORDER
      IF(M .LE. 1) GO TO 9005
      MP = M+1
      JJ = 1
C                                  INITIALIZE WORK VECTOR
      IWK(1) = 1
      DO 75  I = 2,MP
         IWK(I) = IWK(I-1) * 2
   75 CONTINUE
      N4 = IWK(MP-2)
      IF (M .GT. 2) N8 = IWK(MP-3)
      N2 = IWK(MP-1)
      LM = N2
      NN = IWK(MP)+1
      MP = MP-4
C                                  DETERMINE INDICES AND SWITCH A
      J = 2
   80 JK = JJ + N2
      AK2 = A(J)
      A(J) = A(JK)
      A(JK) = AK2
      J = J+1
      IF (JJ .GT. N4) GO TO 85
      JJ = JJ + N4
      GO TO 105
   85 JJ = JJ - N4
      IF (JJ .GT. N8) GO TO 90
      JJ = JJ + N8
      GO TO 105
   90 JJ = JJ - N8
      K = MP
   95 IF (IWK(K) .GE. JJ) GO TO 100
      JJ = JJ - IWK(K)
      K = K - 1
      GO TO 95
  100 JJ = IWK(K) + JJ
  105 IF (JJ .LE. J) GO TO 110
      K = NN - J
      JK = NN - JJ
      AK2 = A(J)
      A(J) = A(JJ)
      A(JJ) = AK2
      AK2 = A(K)
      A(K) = A(JK)
      A(JK) = AK2
  110 J = J + 1
C                                  CYCLE REPEATED UNTIL LIMITING NUMBER
C                                    OF CHANGES IS ACHIEVED
      IF (J .LE. LM) GO TO 80
C
 9005 RETURN
      END

