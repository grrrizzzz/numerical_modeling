      SUBROUTINE SST(T)

c routine for updating ssta; grid is same as atmosphere model, and different
c from ocean dynamics grid.... computations are done dimensionally with
c SSTA expressed in degrees C.  The code is specific to a timestep of 10days
c and other time steps cannot be used correctly without changes to various
c parameters....

C      DIMENSION HBAR(34),TD(30,34),DTZM(34),DTZP(30,34)
      DIMENSION HBAR(30,34),TD(30,34),DTZM(34),DTZP(30,34)
      DIMENSION DT(30,34,7),FIX(34),HZC1(30,34),HDELTA(30,34)
      DIMENSION WM1T(30,34),UV1T(30,34),UV2T(30,34),TT0(30,34)

      COMMON/TIME1/IT,MP,TY

      COMMON/SSTPAR/GAM1,GAM2,TDA1,TDB1,TDA2,TDB2,TLOSS

      COMMON/ZDATA/WM(30,34,2,12),DIVM(30,34,12),Q0O(30,34),
     A SSTM(30,34,12),UO(30,34),VO(30,34),DO(30,34),HTAU(34,30,2),
     B WEM(30,34,12),TO(30,34),UV(30,34,2,12),US(30,34),VS(30,34),
     C WP(30,34),DT1(30,34,7),TT(30,34),UV1(30,34),UV2(30,34),
     D WM1(30,34),UAT(30,34),VAT(30,34),DIVT(30,34),HZC(30,34,12),
     F HLF(30,34,12),TAUMX(1020,1020),TAUMY(1020,1020)
 
      COMMON/ZDAT2/H1(30,34),U1(30,34),V1(30,34)
      COMMON/HMEAN/HM(30,34,12)
      DIMENSION HBAR0(34)

      DATA DTZM/15*.2,.175,.15,.125,.2,.25,.3,.35,.4,.45,.5,9*.55/

      DATA HBAR0/15*1.65,1.7,3*1.75,1.5,1.25,1.,.932,.864,.796,.728,
     + .66,.592,.524,5*.5/

      DATA FIX/15*.15,0.14,3*.13,.2,.3,.45,.48,.5,.55,9*.6/
 
      SAVE
 
1234  DO 6663 I=1, 30
      DO 6663 J=1, 34
      DO 6663 N=1, 7
      DT1(I,J,N)=0.
6663  CONTINUE
 
      UFACTR=5.4
      R=.23
 
c calculate the total temperature and mean currents:
      DO 61 I=5, 26
      DO 62 J=6, 32
      HDELTA(I,J)=HLF(I,J,IT)+TY*(HLF(I,J,MP)-HLF(I,J,IT))
      WM1(I,J)=WEM(I,J,IT)+TY*(WEM(I,J,MP)-WEM(I,J,IT))
      UV1(I,J)=UV(I,J,1,IT)+TY*(UV(I,J,1,MP)-UV(I,J,1,IT))
      UV2(I,J)=UV(I,J,2,IT)+TY*(UV(I,J,2,MP)-UV(I,J,2,IT))
      TT0(I,J)=SSTM(I,J,IT)+TY*(SSTM(I,J,MP)-SSTM(I,J,IT)) ! Other ADV scheme
      TT(I,J)=TO(I,J)+SSTM(I,J,IT)+TY*(SSTM(I,J,MP)-SSTM(I,J,IT))
c      TT(I,J)=SSTM(I,J,IT)+TY*(SSTM(I,J,MP)-SSTM(I,J,IT))  ! Uncoupled
      HBAR(I,J)=(HM(I,J,IT)+TY*(HM(I,J,MP)-HM(I,J,IT)))*.01
62    CONTINUE
61    CONTINUE

c*************************************************
c Linearization of the Tsub function
c
      DO 20 I=6, 25
      DO 25 J=6, 32
      HBARX = HBAR0(J) + HBAR(I,J)
      TD(I,J)=FIX(J)*12.0*(0.5/HBARX)*
     &        (TANH((H1(I,J)-15.0)/35.0)-TANH(-15.0/35.0))
********* AXEL TSUB MODIFICATION
      DTZP(I,J)=.1*(TO(I,J)-TD(I,J))

25    CONTINUE
20    CONTINUE

c      DO 20 I=6, 25
c      DO 25 J=6, 32
c      TD(I,J)=FIX(J)*12.0*(0.5/(HBAR(J)+.01*HDELTA(I,J)))*(TANH(
c     &      (H1(I,J)+HDELTA(I,J)-15.0)/35.0)-TANH(-15.0/35.0))
c
c      DTZP(I,J)=.1*(TO(I,J)-TD(I,J))
c
c25    CONTINUE
c20    CONTINUE
c
c*************************************************

c  compute the surface current anomalies.....constants are such as to give
c  results in cm/sec.  ufactr assumes MLD and upper layer depths of 50m and
c  150m, and the factor of 250. multiplying u1,v1 is because the currents 
c  from mloop.f are coming non-dimensional and they have to be coverted
c  to cm/sec for the coupling.

22    DO 30 I=5, 26
      I1=31-I
      Y=-FLOAT(I-1)*.2+2.9
      A1=UFACTR/(Y*Y+R*R)
      DO 35 J=6, 32
      US(I,J)=U1(I,J)*250.+A1*(R*HTAU(J,I1,1)+Y*HTAU(J,I1,2))
      VS(I,J)=V1(I,J)*250.+A1*(R*HTAU(J,I1,2)-Y*HTAU(J,I1,1))
35    CONTINUE
30    CONTINUE

c next, calculate the anomaly in upwelling velocity
c computation assumes currents are calculated on the 2 X 5.635 degree grid.
 
      DO 42 I=1, 30
      DO 42 J=1, 34
      WP(I,J)=0.
42    CONTINUE
      DX=.5625*2.
      DY=.2*2.
      DO 43 I=6, 25
      DO 44 J=7, 31
      WP(I,J)=.045*(US(I,J+1)-US(I,J-1))/DX+.045*(VS(I-1,J)-
     A VS(I+1,J))/DY
44    CONTINUE
      WP(I,32)=.045*(-US(I,31))/DX+.045*(VS(I-1,32)-VS(I+1,32))/DY
      WP(I,6)=WP(I,7)
43    CONTINUE
 
c calculate the change in sst, using variable sub-timestep to insure
c numerical stability

      CALL GETN(NSST)
c      WRITE(6,3554) NSST
3554  FORMAT(1X, 'NUMBER OF SST LOOPS IS: ', I2)

      DO 50 N=1, NSST
      DO 51 I=6, 25
      DO 52 J=6, 32
      IF(J.EQ.32.OR.J.EQ.6) GO TO 45

c zonal advection
c      DT(I,J,1)=.008*ADV(1,US,TT,I,J)
c      DT(I,J,1)=DT(I,J,1)+.08*ADV(1,UV1,TO,I,J)
c      DT(I,J,2)=0.
c      GO TO 46
c45    DT(I,J,1)=0.
c      DT(I,J,2)=0.

      DT(I,J,1)=.008*ADV(1,UV1T,TT,I,J)        ! Other ADV scheme
      DT(I,J,1)=DT(I,J,1)-.08*ADV(1,UV1,TT0,I,J)
      DT(I,J,2)=0.
      GO TO 46                           
45    DT(I,J,1)=0.
      DT(I,J,2)=0.

c meridional advection
c46    DT(I,J,3)=.008*ADV(2,VS,TT,I,J)
c      DT(I,J,3)=DT(I,J,3)+.08*ADV(2,UV2,TO,I,J)
c      DT(I,J,4)=0.
c      A1=WM1(I,J)
c      A2=WP(I,J)

46    DT(I,J,3)=.008*ADV(2,UV2T,TT,I,J)        ! Other ADV scheme
      DT(I,J,3)=DT(I,J,3)-.08*ADV(2,UV2,TT0,I,J)
      DT(I,J,4)=0.
      A1=WM1(I,J)
      A2=WP(I,J)

c upwelling
      DT(I,J,5)=-.864*HF(A1)*DTZP(I,J)*GAM1
      DT(I,J,6)=-.864*GF(A1,A2)*(DTZM(J)+DTZP(I,J))*GAM2

c surface heat flux
      DT(I,J,7)=(TLOSS-1.)*TO(I,J)*4.
52    CONTINUE
51    CONTINUE

c add up the different tendencies to update ssta
      DO 53 J=6, 32
      DO 54 I=6, 25
      TO(I,J)=TO(I,J)+(DT(I,J,1)+DT(I,J,3)+DT(I,J,5)+DT(I,J,6)+
     + DT(I,J,7)+DT(I,J,4))/FLOAT(NSST)
54    CONTINUE
53    CONTINUE

c update total sst and do not allow anomalies that cause total sst to
c exceed 30 degrees C.
      DO 70 J=6, 32
      DO 71 I=5, 26
c      TT(I,J)=SSTM(I,J,IT)+TY*(SSTM(I,J,MP)-SSTM(I,J,IT))  ! Uncoupled
      TT(I,J)=TO(I,J)+SSTM(I,J,IT)+TY*(SSTM(I,J,MP)-SSTM(I,J,IT))
      IF(TT(I,J) .LE. 30.) GO TO 71
      TO(I,J)=TO(I,J)+30.-TT(I,J)
      TT(I,J)=30.
71    CONTINUE

c update the vertical temp. gradient
      DO 72 I=6, 25
      DTZP(I,J)=.1*(TO(I,J)-TD(I,J))

c save the temperature tendency terms
      DO 73 M=1, 7
      DT1(I,J,M)=DT1(I,J,M)+DT(I,J,M)/FLOAT(NSST)
73    CONTINUE
72    CONTINUE
70    CONTINUE
50    CONTINUE
 
c  save total currents in uv1,uv2,wm1
      DO 81 J=6, 32
      DO 82 I=6, 25
      WM1(I,J)=WM1(I,J)+WP(I,J)
82    CONTINUE
      DO 83 I=5, 26
      UV1T(I,J)=10.*UV1(I,J)+US(I,J)
      UV2T(I,J)=10.*UV2(I,J)+VS(I,J)
83    CONTINUE
81    CONTINUE

C     NOW GET THE U1,V1,H1 FIELDS FOR THE NEXT SST ITERATION
 
215   CALL ZAVG(H1,U1,V1)

      RETURN
      END


      SUBROUTINE GETN(NSST)

c routine to determine number of subintervals to divide 10 days into in
c order to assure numerical stability by CFL type criterion

      INCLUDE 'zeq.common'

      COMMON/ZDATA/WM(30,34,2,12),DIVM(30,34,12),Q0O(30,34),
     A SSTM(30,34,12),UO(30,34),VO(30,34),DO(30,34),HTAU(34,30,2),
     B WEM(30,34,12),TO(30,34),UV(30,34,2,12),US(30,34),VS(30,34),
     C WP(30,34),DT1(30,34,7),TT(30,34),UV1(30,34),UV2(30,34),
     D WM1(30,34),UAT(30,34),VAT(30,34),DIVT(30,34),HZC(30,34,12),
     F HLF(30,34,12),TAUMX(1020,1020),TAUMY(1020,1020)

      COMMON/SGRID/XWRS,DXRS,YNRS,DYRS,NYPS
      SAVE

      DX1=DXRS*1.11E7
      DY1=DYRS*1.11E7
      DTT=DTD*2.628E6
      NSST=2

      DO 10 I=6, 25
      DO 11 J=6, 32
      AAA=ABS(US(I,J))*DTT/DX1
      BBB=ABS(VS(I,J))*DTT/DY1
      IAA=1+AAA
      IBB=1+BBB
      NSST=MAX0(NSST,IAA)
      NSST=MAX0(NSST,IBB)
11    CONTINUE
10    CONTINUE

      print*,nsst

      RETURN
      END


      SUBROUTINE ZAVG(H1,U1,V1)

c routine to average fields from ocean dynamics grid onto the coarser sst grid
      INCLUDE 'zeq.common'
      DIMENSION H1(30,34),U1(30,34),V1(30,34)
      DO 10 I=6, 25

      I1=112-4*I            ! Luis      

      IS=I1-2
      IE=I1+2
      DO 20 J=6, 32
      JJ=J-4                   ! Luis 
      AA=0.
      AA1=0.
      AA2=0.
      DO 21 II=IS, IE
      AA=AA+HB(II,JJ)
      AA1=AA1+UB(II,JJ)
      AA2=AA2+V(II,JJ)
21    CONTINUE
      H1(I,J)=AA/5.0          ! Luis (Div by 15)
      U1(I,J)=AA1/5.0
      V1(I,J)=AA2/5.0
20    CONTINUE
      H1(I,33)=H1(I,32)
      H1(I,34)=H1(I,32)
      U1(I,33)=0.
      U1(I,34)=0.
      V1(I,33)=0.
      V1(I,34)=0.
      DO 15 J=1, 5
      H1(I,J)=H1(I,6)
15    CONTINUE
10    CONTINUE
      DO 12 J=1, 34
      DO 13 I=1, 5
      H1(I,J)=0.
      H1(31-I,J)=0.
      U1(I,J)=0.
      U1(31-I,J)=0.
      V1(I,J)=0.
      V1(31-I,J)=0.
13    CONTINUE
12    CONTINUE
 
      RETURN
      END

c*************************************************
c Linearization of the function HF(X) and GF(X)
c
      REAL FUNCTION HF(X)
      HF=0.5*X*(1.0+tanh(1.0e+2*X))
      RETURN
      END

      REAL FUNCTION GF(X,Y)
      Z=X+Y
      GF=0.5*(Z*(1.0+tanh(1.0e+2*Z))-X*(1.0+tanh(1.0e+2*X)))
      RETURN
      END

c************************************************* 
 
      REAL FUNCTION ADV(ID1,U,T,I,J)

c routine for evaluating advection by modified upwind differencing scheme
      DIMENSION U(30,34),T(30,34)
      DX=2.*.5625
      DY=2.*.2
      IF(ID1.EQ.2) GO TO 11
      IF(U(I,J).GE.0.) ADV=(U(I,J-1)+U(I,J))*(T(I,J-1)-T(I,J))/DX
      IF(U(I,J).LT.0.) ADV=(U(I,J)+U(I,J+1))*(T(I,J)-T(I,J+1))/DX
      RETURN
11    IF(U(I,J).GE.0.) ADV=(U(I+1,J)+U(I,J))*(T(I+1,J)-T(I,J))/DY
      IF(U(I,J).LT.0.) ADV=(U(I-1,J)+U(I,J))*(T(I,J)-T(I-1,J))/DY
      RETURN
      END


