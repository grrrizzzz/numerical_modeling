       SUBROUTINE ATMOS(T)
 
c atmosphere model; steady linear dynamics with nonlinear convergence feedback 
c as described in Zebiak,1986 MWR, and Zebiak and Cane, 1987 MWR.

      COMMON/ZDATA/WM(30,34,2,12),DIVM(30,34,12),Q0O(30,34),
     A SSTM(30,34,12),UO(30,34),VO(30,34),DO(30,34),HTAU(34,30,2),
     B WEM(30,34,12),TO(30,34),UV(30,34,2,12),US(30,34),VS(30,34),
     C WP(30,34),DT1(30,34,7),TT(30,34),UV1(30,34),UV2(30,34),
     D WM1(30,34),UAT(30,34),VAT(30,34),DIVT(30,34),HZC(30,34,12),
     F HLF(30,34,12)

      COMMON/ATMPAR/ALPHA,EPS,BETA,CTOL,IMIN,IMAX,ISTEP
      LOGICAL ISTEP

      COMMON/TIME1/IT,MP,TY

      DIMENSION DIF(16,34),Q0(30,34),RDIV(30,34)
      REAL QM(30,34),QF(30,34)
      COMPLEX Q(80,64),VEC(64),A1,RHS(80),V1(80),C2,AQ,EPS1,A2,B2
      COMPLEX D2,U1(80),U(30,64),VEC1(64),E1(80),E2(80),DIV(30,64)
      INTEGER IWK(534)
      LOGICAL IPRINT
      REAL NINO3
      
      SAVE

      DATA DIF/544*0./
      DATA E1, E2/80*(0.,0.), 80*(0.,0.)/
      DATA Q/5120*(0.,0.)/
      DATA UBAR/0./

c some basic parameter specifications...
      IPRINT=.FALSE.
      NMAX=80
      MMAX=64
      PI=3.1415926535897932
      H1=1.0
      DELTAY=0.2
      YT=7.9
      NMAX1=NMAX-1
 
c compute the paramterized atmospheric heating that is related to ssta.

      CALL ZQGEN(Q0)

c determine whether full anomaly field is being recalculated, or whether
c the anomaly field is only being incremented, based on the NINO3 criterion.

401   ISTEP=.TRUE.
      NINO3=0.
      DO 701 I=13, 18
      DO 701 J=20, 30
      NINO3=NINO3+TO(I,J)/66.
701   CONTINUE
      IF(ABS(NINO3) .LE. 0.1) ISTEP=.FALSE.
 
c if reinitializing, set all the "previous" time arrays to zero.
c      IF(ISTEP) GO TO 412
4111  DO 411 I=1, 30
      DO 411 J=1, 34
      Q0O(I,J)=0.
      UO(I,J)=0.
      VO(I,J)=0.
      DO(I,J)=0.
      QF(I,J)=0.
411   CONTINUE
412   CONTINUE
 
      DO 304 I=6, 25
      DO 305 J=6, 32
      Q0O(I,J)=Q0(I,J)-Q0O(I,J)
      QM(I,J)=Q0O(I,J)*ALPHA
305   CONTINUE
304   CONTINUE

c set the convergence feedback loop counter to zero
      IC=0
      GO TO 380

c top of convergence feedback loop......................................
c here we compute the convergence feedback part of the heating at the
c current iteration

c*************************************************
c Linearization of the convergence function
c
310   DO 21 I=6, 25
      Y=-FLOAT(I+25)*0.2+8.1
      DO 22 J=6, 32
      DD=1.1*(DIVM(I,J,IT)+TY*(DIVM(I,J,MP)-DIVM(I,J,IT)))+DO(I,J)
      DIVT(I,J)=DD*0.9
      QM(I,J)=0.5*BETA*((DD+RDIV(I,J))*(-1.0+
     &     tanh(1.0e+3*(DD+RDIV(I,J))))-
     &     DD*(-1.0+tanh(1.0e+3*DD)))+Q0O(I,J)*ALPHA
22    CONTINUE
21    CONTINUE

c*************************************************

c Now the heating is specified; set up the atmosphere dynamics calculation
380   DO 31 I=1, 80
      DO 31 J=1, 64
      Q(I,J)=(0.,0.)
31    CONTINUE
      DO 231 I=6, 25
      DO 231 J=6, 32
      Q(I+25,J)=QM(I,J)*(1.,0.)
231   CONTINUE

c Fourier transform the heating......
      DO 40 I=31, 50
      DO 50 J=1, MMAX
      VEC(J)=Q(I,J)
50    CONTINUE
      CALL FFT2C(VEC,6,IWK)
      DO 60 J=1, MMAX
      Q(I,J)=VEC(J)
60    CONTINUE
40    CONTINUE

c now solve the tridiagonal system of equations derived by finite differencing
c the equation for the Fourier transform of v, which is separable for each
c wavenumber k.

      DO 70 K=1, 33
      SCALE=-FLOAT(K-1)*PI/18.
      EPS1=EPS*(1.,0.)+(0.,1.)*SCALE*UBAR
      A1=DELTAY*DELTAY*((0.,1.)*SCALE/2./EPS1-EPS1*EPS1/H1-
     1 SCALE*SCALE)
      B=DELTAY/2.
      AQ=(0.,1.)*SCALE*DELTAY*DELTAY/2./EPS1
      DO 80 I=1, 79
      E1(I)=(1.,0.)
      E2(I)=(1.,0.)
      Y=-FLOAT(I)*.2+8.1
      U1(I)=A1-(2.,0.)-B*B*Y*Y*(1.,0.)/H1
80    CONTINUE
      U1(80)=A1-(2.,0.)-B*B*YT*YT*(1.,0.)/H1
      RHS(1)=(AQ*YT*Q(1,K)+B*Q(2,K))/H1
      RHS(80)=(-B*Q(79,K)+AQ*(-YT)*Q(80,K))/H1
      DO 90 I=2, 79
      Y=-FLOAT(I)*0.2+8.1
      RHS(I)=(-B*(Q(I-1,K)-Q(I+1,K))+AQ*Y*Q(I,K))/H1
90    CONTINUE

c the actual tridiagonal solver
      CALL TRID(U1, E1, E2, RHS, NMAX1, NMAX, V1)

c now compute the u and divergence transforms from v 
      A2=EPS1*EPS1+SCALE*SCALE*H1
      B2=0.5*EPS1
      C2=(0.,1.)*H1*SCALE/2./DELTAY
      D2=(0.,1.)*SCALE

      DO 100 I=25, 56
      Y=-FLOAT(I)*0.2+8.1
      U1(I)=(B2*Y*V1(I)-C2*(V1(I+1)-V1(I-1))+D2*Q(I,K))/A2
100   CONTINUE

      DO 110 I=25, 56
      Q(I,K)=V1(I)
110   CONTINUE
      DO 111 I=26, 55
      U(I-25,K)=U1(I)
111   CONTINUE

70    CONTINUE
      DO 81 K=1, 33
      DO 82 I=26, 55
      SCALE=-FLOAT(K-1)*PI/18.
      DIV(I-25,K)=SCALE*(0.,1.)*U(I-25,K)+(Q(I-1,K)-Q(I+1,K))/2./DELTAY
82    CONTINUE
81    CONTINUE

c now invert divergence transform back to physical space
      DO 73 I=26, 55
      DO 75 K=2, 32
      U(I-25,66-K)=CONJG(U(I-25,K))
      Q(I,(66-K))=CONJG(Q(I,K))
      DIV(I-25,66-K)=CONJG(DIV(I-25,K))
75    CONTINUE
73    CONTINUE
      DO 250 I=1, 30
      DO 251 J=1, MMAX
      VEC(J)=CONJG(DIV(I,J))/FLOAT(MMAX)
251   CONTINUE
      CALL FFT2C(VEC,6,IWK)
      DO 255 J=1, 34
      RDIV(I,J)=REAL(VEC(J))
255   CONTINUE
250   CONTINUE

c check to see if the maximum difference of divergence between present and
c preceding iteration falls below the threshold value

      DIFMX=0.0
      DO 501 I=1, 16
      DO 501 J=6, 32
      DIF(I,J)=ABS(RDIV(I+7,J)-DIF(I,J))
      IF(DIF(I,J) .GT. DIFMX) DIFMX=DIF(I,J)
501   CONTINUE
      IF(DIFMX .LE. CTOL .AND. IC .GE. IMIN) IPRINT=.TRUE.
      DO 502 I=1, 16
      DO 502 J=6, 32
      DIF(I,J)=RDIV(I+7,J)
502   CONTINUE

c if another iteration is required and maximum iteration number not exceeded, 
c loop back.
      IF(.NOT. IPRINT .AND. IC .NE. IMAX) THEN
      IC=IC+1
      GO TO 310
      ENDIF

c otherwise, finish up by inverting wind field transforms
      DO 120 I=26, 55
      DO 130 J=1, MMAX
      VEC(J)=(CONJG(U(I-25,J)))/FLOAT(MMAX)
      VEC1(J)=(CONJG(Q(I,J)))/FLOAT(MMAX)
130   CONTINUE
      CALL FFT2C(VEC,6,IWK)
      DO 132 J=1, 34
      U(I-25,J)=VEC(J)
132   CONTINUE
      CALL FFT2C(VEC1,6,IWK)
      DO 134 J=1, 34
      Q(I,J)=VEC1(J)
134   CONTINUE
120   CONTINUE

      DO 135 J=1, 64
      DO 136 I=26, 55
      DIV(I-25,J)=Q(I,J)
136   CONTINUE
135   CONTINUE

c solution for u,v,divergence,ssta-related heating anomaly,and full htg
c anomaly placed in uo,vo,do,q0o,qf, respectively.
257   DO 2001 I=1, 30
      DO 2002 J=1, 34
      UO(I,J)=REAL(U(I,J))+UO(I,J)
      VO(I,J)=REAL(DIV(I,J))+VO(I,J)
2002  CONTINUE
2001  CONTINUE
      DO 2003 I=6, 25
      DO 2004 J=6, 32
      DO(I,J)=RDIV(I,J)+DO(I,J)
      Q0O(I,J)=Q0(I,J)
      QF(I,J)=QM(I,J)+QF(I,J)
2004  CONTINUE
2003  CONTINUE

      ICWR=IC+1
      IF(IC .GT. 0) WRITE(6,2066) T, ICWR
2066  FORMAT(1X,'TIME=',F6.2,5X,'NO. OF ITERATIONS:',1X, I1)
      IF(.NOT. ISTEP) WRITE(6,2067) T
2067  FORMAT(1X, 'REINITIALIZED AT T=',F7.2)

      RETURN
      END

      SUBROUTINE ZQGEN(Q0)

      DIMENSION Q0(30,34),H1(30,34)

      COMMON/TIME1/IT,MP,TY

      COMMON/ZDATA/WM(30,34,2,12),DIVM(30,34,12),Q0O(30,34),
     A SSTM(30,34,12),UO(30,34),VO(30,34),DO(30,34),HTAU(34,30,2),
     B WEM(30,34,12),TO(30,34),UV(30,34,2,12),US(30,34),VS(30,34),
     C WP(30,34),DT1(30,34,7),TT(30,34),UV1(30,34),UV2(30,34),
     D WM1(30,34),UAT(30,34),VAT(30,34),DIVT(30,34),HZC(30,34,12),
     F HLF(30,34,12)
C
      COMMON/CPAR/FATM,AMPL

      SAVE
 
c first smooth the ssta field....
      DO 2011 I=6, 25
      DO 2011 J=6, 32
      Q0(I,J)=TO(I,J)
2011  CONTINUE
      DO 201 N=1, 3
       DO 202 J=6, 32
      DO 203 I=7, 24
      H1(I,J)=(Q0(I-1,J)+2.*Q0(I,J)+Q0(I+1,J))/4.
203   CONTINUE
      DO 205 I=1, 6
      I1=31-I
      H1(I,J)=H1(7,J)*FLOAT(I-1)*.166
      H1(I1,J)=H1(24,J)*FLOAT(I-1)*.166
205   CONTINUE
      DO 204 I=6, 25
      Q0(I,J)=H1(I,J)
204   CONTINUE
202   CONTINUE
201   CONTINUE
      DO 210 I=6, 25
      DO 211 J=7, 31
      H1(I,J)=(Q0(I,J-1)+Q0(I,J)*2.+Q0(I,J+1))/4.
211   CONTINUE
      H1(I,6)=(TO(I,6)*3.+TO(I,7))/4.
      H1(I,32)=(TO(I,31)+3.*TO(I,32))/4.
      DO 212 J=6, 32
      Q0(I,J)=H1(I,J)
212   CONTINUE
210   CONTINUE
 
C          DO THE Q0 PARAMETERIZATION

5556  DO 5555 I=6, 25
      DO 55 J=6, 32
      A=SSTM(I,J,IT)+TY*(SSTM(I,J,MP)-SSTM(I,J,IT))
      Q0(I,J)=Q0(I,J)*EXP(0.06*(A-29.8))
55    CONTINUE
5555  CONTINUE
 
      RETURN
      END
 
 
      SUBROUTINE STRESS(T,ISTART)

c routine for computing stress for ocean model
C      IF NIC=1, A PSEUDO RANDOM W. PAC. STRESS FORCING IS ADDED 
C     IF NIC=4, AN INITIAL CONDITION WIND IS USED FOR T<TCUT, GOVERNED 
C     BY TAMP1,ITOP1,IBOT1,JLEFT1,JRT1, UNIFORM IN SPACE

      COMMON/TIME1/IT,MP,TY

      COMMON/ZDATA/WM(30,34,2,12),DIVM(30,34,12),Q0O(30,34),
     A SSTM(30,34,12),UO(30,34),VO(30,34),DO(30,34),HTAU(34,30,2),
     B WEM(30,34,12),TO(30,34),UV(30,34,2,12),US(30,34),VS(30,34),
     C WP(30,34),DT1(30,34,7),TT(30,34),UV1(30,34),UV2(30,34),
     D WM1(30,34),UAT(30,34),VAT(30,34),DIVT(30,34),HZC(30,34,12),
     F HLF(30,34,12)

      double precision DSEED
      DIMENSION UT(30,34)
      REAL NINO3

      COMMON/STRPAR/TCUT,TSEED,TAMP1,TAMP2,JLEFT,JRT,ITOP,IBOT,NIC
      COMMON/CPAR/FATM,AMPL

      SAVE

      IF (ISTART .EQ. 0) THEN
      DSEED=TSEED*1.
c      RNORM=GGNQF(DSEED)      
c      RNORM=GGNQF(DSEED)
      PI=3.1415926535897932
      OMEGA=3*PI/2.
c      T0=2.*PI*GGUBFS(DSEED)
c      RNORM=GGNQF(DSEED)
      R=.5
      ENDIF

      DO 5 I=1, 30
      DO 5 J=1, 34
      UT(I,J)=0.
5     CONTINUE
 

C        DEFINE UTEMP HERE

      IF(NIC .EQ. 4 .AND. T .LE. TCUT) THEN
      DO 17 I=ITOP, IBOT
      DO 18 J=JLEFT, JRT
      UT(I,J)=TAMP1
18    CONTINUE
17    CONTINUE
      ELSE IF(NIC .EQ. 1) THEN
      NINO3=0.
      DO 701 I=13, 18
      DO 701 J=20, 30
      NINO3=NINO3+TO(I,J)/66.
701   CONTINUE
c      RNORM=GGNQF(DSEED)
      AMPR=(TAMP1+TAMP2*RNORM)*COS(OMEGA*T+T0)
      IF(NINO3 .GT. 1.) AMPR=AMPR/SQRT(NINO3)
21    DO 15 I=ITOP, IBOT
      Y=-FLOAT(I-1)*.2+2.9
      DO 16 J=JLEFT, JRT
      X=5.625*(J-(JLEFT1+JRT1)/2.)/30.
      UT(I,J)=R*UT(I,J)+AMPR*EXP(-X*X)*EXP(-Y*Y)
16    CONTINUE
15    CONTINUE

      ENDIF

c now compute the stress from winds using bulk formula; large drag coeff.
c Change the coeff from 0.2 to 2.0 in the sensitivity studies. The standard
c coeff in HTAU is 1.0

      cc=FATM
 
11    DO 10 I=6, 25
      I1=31-I
      XLAT=-29.+2.*(I1-1)
      DO 20 J=6, 32
      XLON=101.25+5.625*(J-1)
      A1=UO(I,J)+WM(I,J,1,IT)+TY*(WM(I,J,1,MP)-WM(I,J,1,IT))
      UAT(I,J)=A1
      A2=VO(I,J)+WM(I,J,2,IT)+TY*(WM(I,J,2,MP)-WM(I,J,2,IT))
      VAT(I,J)=A2
      AT=SQRT(A1*A1+A2*A2)
      A3=WM(I,J,1,IT)+TY*(WM(I,J,1,MP)-WM(I,J,1,IT))
      A4=WM(I,J,2,IT)+TY*(WM(I,J,2,MP)-WM(I,J,2,IT))
      AM=SQRT(A3*A3+A4*A4)
C      HTAU(J,I1,1)=(A1*AT-A3*AM)*.0329+UT(I,J)
      HTAU(J,I1,1)=cc*(A1*AT-A3*AM)*.0329
     &      +AMPL*EXP(-((XLON-160.)/20.)**2-(XLAT/10.)**2)
      HTAU(J,I1,2)=cc*(A2*AT-A4*AM)*.0329

20    CONTINUE
10    CONTINUE

c taper the northern and southern regions to avoid spurious wave generation

      DO 30 J=6, 32
      DO 31 K=1, 2
      DO 32 I=2, 5
      I1=31-I
      HTAU(J,I,K)=0.2*FLOAT(I-1)*HTAU(J,6,K)       
      HTAU(J,I1,K)=0.2*FLOAT(I-1)*HTAU(J,25,K)    
32    CONTINUE
31    CONTINUE
30    CONTINUE

      RETURN
      END
 
      SUBROUTINE TRID(D, DP1, DM1, B, N1, N, X)

c tridiagonal solver
      INTEGER N
      COMPLEX D(80), DP1(80), DM1(80), B(80), X(80)
      INTEGER I, N1
      DP1(1)=DP1(1)/D(1)
      DO 170 I=2, N1
      D(I)=D(I)-DM1(I-1)*DP1(I-1)
      DP1(I)=DP1(I)/D(I)
170   CONTINUE
      D(N)=D(N)-DM1(N1)*DP1(N1)
      B(1)=B(1)/D(1)
      DO 180 I=2, N
      B(I)=(B(I)-DM1(I-1)*B(I-1))/D(I)
180   CONTINUE
      X(N)=B(N)
      DO 190 I=1, N1
      X(N-I)=B(N-I)-DP1(N-I)*X(N+1-I)
190   CONTINUE
      RETURN
      END


