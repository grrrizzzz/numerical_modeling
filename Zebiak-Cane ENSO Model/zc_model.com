#!/bin/csh -x


foreach GAM (0.75)
set GAM1=${GAM}
set GAM2=${GAM}
#### COUPLING COEFFICIENT
#foreach  MU (0.9 1.0 1.1 1.2 1.3)
set MU=0.95

set FC=/opt/intel_fc_80/bin/ifort
set FC=g77
#set FFLAG0='-r8 -i4'
echo 'hier'
#set FFLAG1='-O2'
#set FFLAG0=''
#set FFLAG1=''

cat>fctest.data<<EOF
Test run of coupled ZC model
zeq9fsu.hst
outhst
NSTART =  0
TFIND  =  120.5
NX     =  30
NY     =  100
HEQUIV =  64.
TZERO  =  0.5
TENDD  =  1200.5
DTD    =  .3333334
XWD    =  124.
XED    =  280.
YSD    =  -28.75
YND    =   28.75
TDECAY =   30.
TPLMIN =  -1000.
NPRINT =  3
NSEG   =  1
YNORTH =  28.75
YSOUTH =  -28.75
XWEST  =  124.
NTAPE  =  0
NREWND =  11
NATM   =  2
IC2    =  0
GAM1   =  ${GAM1}
GAM2   =  ${GAM2}
TDA1   =  28.
TDB1   =  1.25
TDA2   =  -40.
TDB2   =  3.0
TLOSS  =  .98
ALPHA  =  1.6
EPS    =  .3
BETA   =  .75
CTOL   =  .15
IMIN   =  1
IMAX   =  2
ISTEP  = .TURE.
TCUT   =  4.
TSEED  =  4.
TAMP1  =  0.5
TAMP2  =  0.
JLEFT  =  12
JRT    =  18
ITOP   =  13
IBOT   =  18
NIC    =  0
EOF



cat>axmain_cc.f<<EOF
c*****************************************************************
c This is the main program (driver routine) for the coupled atmos*
c ocean model by Zebiak and Cane (MWR 1987). Note time convention*
c in months; 0.5 corresponds to mid January 1960.                *
c*****************************************************************

      INCLUDE 'zeq.common'

      COMMON/ZDATA/WM(30,34,2,12),DIVM(30,34,12),Q0O(30,34),
     A SSTM(30,34,12),UO(30,34),VO(30,34),DO(30,34),HTAU(34,30,2),
     B WEM(30,34,12),TO(30,34),UV(30,34,2,12),US(30,34),VS(30,34),
     C WP(30,34),DT1(30,34,7),TT(30,34),UV1(30,34),UV2(30,34),
     D WM1(30,34),UAT(30,34),VAT(30,34),DIVT(30,34),HZC(30,34,12),
     F HLF(30,34,12),TAUMX(1020,1020),TAUMY(1020,1020)

      COMMON/TIME1/IT,MP,TY
      COMMON/INIT/IC2
      COMMON/TFED/USST(30,34),VSST(30,34),HTAUSST(34,30,2)
      COMMON/DHEAT/Q2(30,34)
      COMMON/CPAR/FATM,AMPL
c
      DIMENSION DUMMY(30,34,14),HH1(30,34),UU1(30,34),VV1(30,34)
c
      REAL NINO3
      REAL*4 SNINO3
      DOUBLE PRECISION NN
      DATA TFIND/1.E20/
c
c Reading the parameters (information) of the couple model
c
      OPEN(5,FILE='fctest.data',status='old')

c Opening the files to save the output data
      NN=1.D0
      FATM0 = 0.7
      IREC = 0
      IREC1 = 0
      NINO3 = 0.
      YNINO = 0.

      open(91,file='nino.gdat'
     *       ,form='formatted')

c      open(92,file='uvh0.gdat'
c     *       ,form='unformatted',
c     *        access='direct',recl=4*30*34,status='unknown')

c      open(92,file='v15.gdat'
c     *       ,form='unformatted',
c     *        access='direct',recl=4*30*34,status='unknown')
c
c Read case title information and write to standard output
c First line in the FILE 5.
c
      READ(5,450) HEADER 
450   FORMAT(80A1)                                                      
      WRITE(6,475) HEADER
475   FORMAT('  CASE TITLE: ',80A1/)
c 
c First, open necessary files to run the coupled model.
c
      CALL OPENFL
c 
c Read control parameters from input data file on unit 5.
c
      NATM=2
      READ(5,40,END=999) NSTART, TFIND
      IF(NSTART.LT.0) GO TO 999
      WRITE(6,50) NSTART, TFIND
40    FORMAT ( 8X, I15 / 8X, F15.0 )
50    FORMAT ( / I15, ' NSTART ' / E15.7, ' TFIND ' / )
c
c If history input is being used for startup, call rdhist.
c     
c      IF(NSTART.GT.0) CALL RDHIST(0,TFIND)
      NATMR=NATM
c
c Now pick up (possibly) new ocean model parameters unless nstart=1.
c
      IF(NSTART.NE.1) CALL SETUP
c 
c If nstart=1 or nstart=2, time is gotten from the history file,
c otherwise (nstart=3) time is set to tzero from input data file.
c
      IF(NSTART .EQ. 3 .OR. NSTART .EQ. 0) THEN
      TD=TZERO
      ELSE
      TZERO=TD
      ENDIF

      NT=0
      IT1=TD+120.501
      IT=MOD(IT1, 12)
      MP=IT+1
      IF(IT.EQ.0) IT=12
      TY=TD+120.501-IT1

c fixed July condition
c      IT=MOD(127.001, 12)
c      MP=IT+1
c      IF(IT.EQ.0) IT=12
c      TY=6.5+120.501-127

c Read in ssta model and atmosphere model constants, etc...

       CALL SETUP2 
c Read in mean field data.  This data was open in the OPENFL call.

       CALL INITDAT

c Do initialization of ocean dynamics model constants, etc...

       CALL CONSTC

c Generating Hermite Functions

       CALL HERMITE

c Setting matrix for ocean model

       CALL MATRIX

       ISTART=0
       
       IF(NSTART .GT. 0) THEN
       ISTART=1
       GO TO 400
       ELSE
       GO TO 350
       ENDIF

c***********************************************************************
c this is the top of main loop for all but the first cycle

c Update the ssta field. Using HTAU, we compute the SST anonaly and
c and total SST (arrays T0 and TT respectively) fields.

300   CALL SST(TD)



*      write(*,*) 'halloo',TD
*      CALL ATMOS(TD)
*      CALL STRESS(TD,ISTART)
c Given the temperature I compute the stress (using the matrixes)

      CALL TAUMATRIX

c Insert the stress forcing in ocean model..............................

      CALL CFORCE
      CALL XMYM
 
c update the ocean dynamics.............................................

350   CALL WAVES
      CALL CURRENTS
C
c      AMPL = 0.7*2.*(GGUBFS(NN)-0.5)
c      AMPL = 0.0007*2.*(GGUBFS(NN)-0.5)
************* FOR WIND-STRESS NOISE
*       AMPL = 0.2*2.*(GGUBFS(NN)-0.5)
C      AMPL = 0.00001*2.*(GGUBFS(NN)-0.5)
      AMPL = 0.0
      NN = NN + 1.D0
C      FATM0 = FATM0 - 0.6/(250.*36.)
C      FATM = 0.837*FATM0
C      FATM = 0.805*FATM0
      FATM = ${MU}
C
400   CONTINUE

c///////////////////////////////////////////////////////////////////////
c  here one can add any desired output write statements.................
c///////////////////////////////////////////////////////////////////////

c write history output if conditions are met............................
c      IF(NTAPE.NE.0.AND.MOD(NT,NTAPE).EQ.0.AND.NT.NE.0.AND.
c     * TD .GE. TPLMIN) CALL WRHIST

c write the results
c      IF(NTAPE.NE.0.AND.MOD(NT,NTAPE).EQ.0.AND.NT.NE.0.AND.
c     * TD .GE. TPLMIN) THEN
c
      XNINO=0.
      DO 701 I=13, 18
      DO 701 J=20, 30
      XNINO=XNINO+TO(I,J)/66.
701   CONTINUE
      MON1=TD
      IF(IREC.EQ.0) MON2=TZERO
c
      IF(MON1 .NE. MON2) THEN
      IREC=IREC+1
      NINO3=NINO3/YNINO
c      SNINO3=SNGL(NINO3)
      SNINO3=(NINO3)
c
      DO 705 I=1,30
      DO 705 J=1,34
      DO 705 NP=1,14
      DUMMY(I,J,NP)=DUMMY(I,J,NP)/YNINO
705   CONTINUE
c
      IF(NINO3 .LT. XMINNINO3) XMINNINO3=NINO3
      IF(NINO3 .GT. XMAXNINO3) XMAXNINO3=NINO3
c      PRINT*, 'MAX= ',XMAXNINO3, '  MIN= ', XMINNINO3
c
*      WRITE(91,REC=IREC) SNINO3
      WRITE(91,910) NINO3
910   format(12F7.2)

c      WRITE(92,REC=5*IREC-4) ((DUMMY(I,J,1),J=1,34),I=30,1,-1)
c      WRITE(92,REC=5*IREC-3) ((DUMMY(I,J,11),J=1,34),I=30,1,-1)
c      WRITE(92,REC=5*IREC-2) ((DUMMY(I,J,12),J=1,34),I=30,1,-1)
c      WRITE(92,REC=5*IREC-1) ((DUMMY(I,J,13),J=1,34),I=30,1,-1)
c      WRITE(92,REC=5*IREC-0) ((DUMMY(I,J,14),J=1,34),I=30,1,-1)

c      WRITE(92,REC=IREC) ((DUMMY(I,J,11),J=1,34),I=30,1,-1)

      NINO3=0.
      YNINO=0.
      DO 706 I=1,30
      DO 706 J=1,34
      DO 706 NP=1,14
      DUMMY(I,J,NP)=0.
706   CONTINUE
      MON2=MON1
      ENDIF
c
      NINO3=NINO3+XNINO
c
      CALL ZAVG(HH1,UU1,VV1)
c
      DO 702 I=1,30
      DO 702 J=1,34
      DUMMY(I,J,1)=DUMMY(I,J,1)+TO(I,J)
      DUMMY(I,J,2)=DUMMY(I,J,2)+HTAU(J,I,1)
      DUMMY(I,J,3)=DUMMY(I,J,3)+HTAU(J,I,2)
      DUMMY(I,J,4)=DUMMY(I,J,4)+DT1(I,J,1)-DT1(I,J,2)
      DUMMY(I,J,5)=DUMMY(I,J,5)+DT1(I,J,2)
      DUMMY(I,J,6)=DUMMY(I,J,6)+DT1(I,J,3)-DT1(I,J,4)
      DUMMY(I,J,7)=DUMMY(I,J,7)+DT1(I,J,4)
      DUMMY(I,J,8)=DUMMY(I,J,8)+DT1(I,J,5)
      DUMMY(I,J,9)=DUMMY(I,J,9)+vv1(I,J)
      DUMMY(I,J,10)=DUMMY(I,J,10)+uu1(I,J)
      DUMMY(I,J,11)=DUMMY(I,J,11)+HH1(I,J)
      DUMMY(I,J,12)=DUMMY(I,J,12)+US(I,J)
      DUMMY(I,J,13)=DUMMY(I,J,13)+VS(I,J)
      DUMMY(I,J,14)=DUMMY(I,J,14)+WP(I,J)

702   CONTINUE
      
      YNINO=YNINO+1.

c 
c check to see if it's time to stop
c
      IF(NT.GE.NTIMES) GO TO 900

c if not, loop again
      NT=NT+1
      TD=TZERO+FLOAT(NT)*DTD
      IT1=TD+120.501
      IT=MOD(IT1, 12)
      MP=IT+1
      IF(IT.EQ.0) IT=12
      TY=TD+120.501-IT1

c fixed july condition
c      IT=MOD(127.001, 12)
c      MP=IT+1
c      IF(IT.EQ.0) IT=12
c      TY=6.5+120.501-127
      GO TO 300

c time to finish up

c900   CALL WRHIST

 900  CONTINUE
 
       WRITE(6,998)
998    FORMAT(' **************** NORMAL END-OF-JOB ZEQ *************')
 
999   CONTINUE

       STOP
       END
C
      REAL FUNCTION GGUBFS (DSEED)
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   D2P31M,D2P31
C                                  D2P31M=(2**31) - 1
C                                  D2P31 =(2**31)(OR AN ADJUSTED VALUE)
      DATA               D2P31M/2147483647.D0/
      DATA               D2P31 /2147483648.D0/
C                                  FIRST EXECUTABLE STATEMENT
      DSEED = DMOD(16807.D0*DSEED,D2P31M)
      GGUBFS = DSEED / D2P31
      RETURN
      END
 

EOF



rm -f axopenfl_FOM.f
cat>axopenfl_FOM.f<<EOF
c********************************************************
c This routine opens the files to run the model.        *
c                                                       *
c rcwindmn.data: Climatological monthly mean surface    *
c                winds.                                 *
c                                                       *
c rcdivmn.data: Climatological monthly mean surface     *
c               wind divergence.                        *
c                                                       *
c rcsstmn.data: Climatological monthly mean sea surface *
c               temperature.                            *
c                                                       *
c wem.zeb: Simulated climatological monthly mean oceanic* 
c          upwelling.                                   *
c                                                       *
c uv.zeb: Simulated climatological monthly mean surface * 
c          horizontal currents.                         *
c                                                       *
c********************************************************

      SUBROUTINE OPENFL

      CHARACTER*45 FN61,FN62
 
      READ(5,111) FN61
111     FORMAT(A45)
      OPEN(UNIT=61, FILE=FN61, FORM='UNFORMATTED')
      REWIND 61
      READ(5,112) FN62
112     FORMAT(A45)
      OPEN(UNIT=62, FILE=FN62, FORM='UNFORMATTED')
C
      OPEN(UNIT=82,FILE='rcwindmn.data_FOM', STATUS='OLD')
      OPEN(UNIT=83,FILE='rcdivmn.data_FOM', STATUS='OLD')
      OPEN(UNIT=84,FILE='rcsstmn.data_FOM', STATUS='OLD')
      OPEN(UNIT=87,FILE='wem.zeb_FOM', STATUS='OLD')
      OPEN(UNIT=89,FILE='uv.zeb_FOM', STATUS='OLD')
      open(53,file='MatrixTx.data',status='old')
      open(54,file='MatrixTy.data',status='old')
c
      RETURN
      END

EOF

rm -f axtaumatrix.f
cat>axtaumatrix.f<<EOF
c*****************************************************************
c This is the main program (driver routine) for the coupled atmos*
c ocean model by Zebiak and Cane (MWR 1987). Note time convention*
c in months; 0.5 corresponds to mid January 1960.                *
c*****************************************************************

      SUBROUTINE TAUMATRIX

      COMMON/ZDATA/WM(30,34,2,12),DIVM(30,34,12),Q0O(30,34),
     A SSTM(30,34,12),UO(30,34),VO(30,34),DO(30,34),HTAU(34,30,2),
     B WEM(30,34,12),TO(30,34),UV(30,34,2,12),US(30,34),VS(30,34),
     C WP(30,34),DT1(30,34,7),TT(30,34),UV1(30,34),UV2(30,34),
     D WM1(30,34),UAT(30,34),VAT(30,34),DIVT(30,34),HZC(30,34,12),
     F HLF(30,34,12),TAUMX(1020,1020),TAUMY(1020,1020)

      COMMON/TIME/IT,MP,TY

      COMMON/INIT/IC2
c
      COMMON/TFED/USST(30,34),VSST(30,34),HTAUSST(34,30,2)
      COMMON/DHEAT/Q2(30,34)
      COMMON/CPAR/ FATM,AMPL
c
      DIMENSION TO1D(1024),VTAUX(1024),VTAUY(1024)  ! Luis
c
c*************************************************
c Putting T0(i,j) as a one-dim array

      k3=1

      do jj=1,34
         do ii=1,30

            TO1D(k3)=TO(ii,jj)
               
            k3=k3+1

         enddo
      enddo
c
c***********************************************
c Multiply the matrices
c
      do k=1,1024
         VTAUX(k)=0.0
         VTAUY(k)=0.0
      enddo


      do k1=1,1020
         do k2=1,1020

            VTAUX(k1)=VTAUX(k1)+TAUMX(k2,k1)*TO1D(k2)/0.1
            VTAUY(k1)=VTAUY(k1)+TAUMY(k2,k1)*TO1D(k2)/0.1
            
         enddo
      enddo

c      cc=0.85                            ! coupling coeff
      cc=FATM                             ! coupling coeff
      k4=1

*       write(*,*) 'hallooooo',cc


      do jj=1,34
         do ii=1,30

            HTAU(jj,ii,1)=cc*VTAUX(k4)
            HTAU(jj,ii,2)=cc*VTAUY(k4)
               
            k4=k4+1

         enddo
      enddo

      RETURN
      END
 

EOF

rm -f axsetup.f
cat>axsetup.f<<EOF
      SUBROUTINE SETUP

c read in the data inputs for the ocean dynamics model.....
 
C******** THE INFO PARAMETERS ARE                              
C                                                                       
C NX     THE NUMBER OF GRID INTERVALS IN THE X DIRECTION (NX+1 = NO. OF 
C        U,H POINTS,NX = NO. OF V POINTS)                               
C NY     NUMBER OF GRID INTERVALS IN THE Y DIRECTION (NY+1 =NO. OF V    
C        POINTS, NY = NO. OF U,H POINTS)                                
C HEQUIV THE EQUIVALENT DEPTH (IN CENTIMETERS)                          
C                                                                       
C        .............................................................. 
C        .                   NOTE                                     . 
C        .ALL TIMES ARE IN MONTHS                                     . 
C        .ALL DISTANCES ARE IN DEGREES (=111 KM)                      . 
C        .HOWEVER, IF HEQUIV =0                                       . 
C        .THEN ALL TIMES AND DISTANCES ARE CONSIDERED NONDIMENSIONAL  . 
C        .............................................................. 
C                                                                       
C TZERO  THE INITIAL TIME                                               
C TENDD  THE FINAL TIME                                                 
C DTD    THE TIMESTEP                                                   
C TDECAY RAYLEIGH FRICTION PARAMETER = SPINDOWN TIME                    
C                                                                       
C XWD    WESTERN MOST LONGITUDE OF THE OCEAN BASIN (A U,H POINT)        
C XED    EASTERN MOST LONGITUDE OF THE BASIN (A U,H POINT)              
C YSD    SOUTHERN MOST LATITUDE OF THE BASIN (A V POINT)                
C YND    NORTHERN MOST LATITUDE OF THE BASIN (A V POINT)                
C........THE NEXT 4 PARAMETERS ARE FOR THE SEGMENTED GEOMETRY           
C NSEG   THE NUMBER OF SEGMENTS                                         
C YSOUTH(K) SOUTHERN LATITUDE OF SEGMENT K (A V POINT)                  
C YNORTH(K) NORTHERN LATITUDE OF SEGMENT K (A V POINT)                  
C XWEST(K)  WESTERN LONGITUDE OF SEGMENT K (A U,H POINT)                
C                                                                       
C TPLMIN FIRST TIME FOR WHICH History OUTPUT MAY BE PRODUCED              
C NPRINT PRINTS TIMESTEP SUMMARY EVERY NPRINT TIMES                     
C                                                                       
C***********************************************************************
C***********************************************************************
      INCLUDE 'zeq.common'
 
      WRITE ( 6, 25 )
   25 FORMAT ( // ' Here are the present values of the "INFO" data' /)
 
      READ (5,40,END=4000) NX,NY,HEQUIV,TZERO,TENDD,DTD,XWD,XED,YSD,YND 
     A             ,TDECAY,TPLMIN,NPRINT,NSEG
      DO 1 N=1, NSEG
      READ (5,41,END=4000) YNORTH(N)
1     CONTINUE
      DO 2 N=1, NSEG
      READ (5,41,END=4000) YSOUTH(N)
2     CONTINUE
      DO 3 N=1,NSEG
      READ (5,41,END=4000) XWEST(N)
3     CONTINUE

      READ (5,42,END=4000) NTAPE,NREWND,NATM

40    FORMAT (
     A 8X, I15 /      8X, I15 /
     F 8X, F15.0 /    8X, F15.0 /
     F 8X, F15.0 /    8X, F15.0 /
     F 8X, F15.0 /    8X, F15.0 /
     F 8X, F15.0 /    8X, F15.0 /
     F 8X, F15.0 /    8X, F15.0 /
     A 8X, I15 /      8X, I15   )

41    FORMAT(8X, F15.0)

42    FORMAT(
     F 8X, I15   /    8X, I15   /
     A 8X, I15   )

      WRITE (6,400)  
C
      WRITE ( 6, 50 )NX,NY,HEQUIV,TZERO,TENDD,DTD,XWD,XED,YSD,YND       
     A             ,TDECAY,TPLMIN,NPRINT,NSEG

      WRITE(6,51) YNORTH(1)
      IF(NSEG .GT. 1) THEN
      DO 4 N=2, NSEG
      WRITE(6,52) YNORTH(N)
4     CONTINUE
      ENDIF

      WRITE(6,53) YSOUTH(1)
      IF(NSEG .GT. 1) THEN
      DO 5 N=2, NSEG
      WRITE(6,52) YSOUTH(N)
5     CONTINUE
      ENDIF

      WRITE(6,54) XWEST(1)
      IF(NSEG .GT. 1) THEN
      DO 6 N=2, NSEG
      WRITE(6,52) XWEST(N)
6     CONTINUE
      ENDIF

      WRITE(6,55) NTAPE,NREWND,NATM

   50 FORMAT ( ' NX     =',I15,   5X, ' NY     =',I15   /
     3         ' HEQUIV =',E15.7, 5X, ' TZERO  =',E15.7 /
     5         ' TENDD  =',E15.7, 5X, ' DTD    =',E15.7 /
     7         ' XWD    =',E15.7, 5X, ' XED    =',E15.7 /
     9         ' YSD    =',E15.7, 5X, ' YND    =',E15.7 /
     B         ' TDECAY =',E15.7, 5X, ' TPLMIN =',E15.7 /
     F         ' NPRINT =',I15,   5X, ' NSEG   =',I15  // )

51    FORMAT( ' YNORTH =',E15.7)
52    FORMAT( '         ',E15.7)
53    FORMAT( ' YSOUTH =',E15.7)
54    FORMAT( ' XWEST  =',E15.7)

55    FORMAT(// ' NTAPE  =',I15,  5X, ' NREWND =',I15 /
     P          ' NATM   =',I15  / )
 
  400 FORMAT('*************************  START OF CASE DOCUMENTATION',
     *    '**************************',/)                               
 
      RETURN
 
 4000 CONTINUE
      WRITE ( 6, 50 )
      WRITE ( 6, 5000 )
 5000 FORMAT ( ' END OF FILE, UNIT 5 ' )
      STOP
      END

EOF

rm -f axconstc.f
cat>axconstc.f<<EOF
      SUBROUTINE CONSTC
 
c compute commonly used constants to be saved

      INCLUDE 'zeq.common'

      BETA = 2.28E-11
      SQRT2=SQRT(2.e0)
      PI =ACOS(-1.e0)

      rho=1026.0
      hterm=150.0

c eleq is in degrees latitude, teq is in days

      CWAVE=SQRT(9.81E-2*HEQUIV)
      ELEQ= SQRT(CWAVE/BETA)/111.E3
      TEQ=1.E0/ SQRT(CWAVE*BETA)/86400.
      TEQM=TEQ*12.E0/365.25

c scale the wind stress forcing so that H will come out in meters
c it is assumed the wind stress is specified in dy/cm/cm

      modes=31
      ndim=modes*(nx-1)
      stress0=rho*hterm*cwave*sqrt(cwave*beta)
      WSCALE=(0.1/stress0)
      WRITE(6,7) HEQUIV,CWAVE,ELEQ,TEQ,TEQM,WSCALE
    7 FORMAT(//' HEQUIV =',F4.0,' CM.  CWAVE =',F5.2,' M/S.  LEQ =',
     A         F5.2,' DEGREES.  TEQ =',F5.2,' DAYS =',F7.4,' MONTHS.'/
     B         ' WSCALE =',F5.2,' M/(CM/SEC)**2'/)
 
      XE=(XED-XWD)/ELEQ
      YS=YSD/ELEQ
      YN=YND/ELEQ
      NTIMES=(TENDD-TZERO)/DTD+0.1
      DX = XE/FLOAT(NX)
      DY = (YN-YS)/FLOAT(NY)
      NXP = NX+1
      NYP = NY+1
      DT=DTD/TEQM
c      RFRIC=TEQM/TDECAY
      RFRIC=0.0
      DXD=DX*ELEQ
      DYD=DY*ELEQ
 
      WRITE(6,8) XWD,XED,YSD,YND,DXD,DYD,DTD,TDECAY,XE,YS,YN,DX,DY,DT
     A     ,RFRIC
    8 FORMAT(25X,'XW        XE        YS        YN        DX       DY'
     A     ,'       DT       TDECAY'
     A     ,/' DIMENSIONAL      ',4F10.2,3F10.5,F10.2
     A     ,/' NONDIMENSIONAL   ','       0.0',3F10.2,4F10.5/)
      WRITE(6,9) NX,NY,TZERO,TENDD,NTIMES
    9 FORMAT(/'  NX =',I5,'  NY =',I5,'  TZERO = ',F7.2,
     A     '  TENDD = ',F10.2,'  NTIMES = ',I7/)
C**********************************************************************
 
      DO J=1,NY
         Y(J)=YS+(J-0.5)*DY
      ENDDO
 
      DO 215 I=1, NXP
         DO 215 J=1, NYP
            XXT(J,I)=0.E0
            YYT(J,I)=0.E0
 215     CONTINUE
            
         RETURN
         END

EOF

rm -f axwaves-1.f
cat>axwaves-1.f<<EOF
c************************************************************
c Ocean model based on Battisti's paper (JAS October 1988)  *
c                                                           *
c************************************************************ 

      subroutine waves      

      include 'zeq.common'
      
      print*,ntimes,nt
c
c Getting the kelvin (q0) and Roosby (qm) waves.  Using the
c VectorX subroutine we save the waves at each time step. 
c
      call VectorQ
      call RHS
      call tridiagonal
      call VectorQ
      call VectorX

      return
      end
c
c******** End of the main program**************
c
c***********************************************
c Solving the Right Hand Side of the Equations *
c***********************************************

      subroutine RHS

      include 'zeq.common'

c      implicit none
      integer i,m
      real cs
      INTEGER NX,NY,NXP,NYP,IALPHA,ISTART,NT,NTIMES,NPRINT
      INTEGER NSEG,NTAPE,NREWND,NSTART,NATM,NATMR,MODES,NDIM
      REAL DT,DX,DY,TPLMIN,HK,PI,SQRT2,TD,TZERO,XE,YN,YS
      REAL OMEGA,ENDT,RALPHA,HEQUIV,XWD,XED,YSD,YND,CWAVE
      REAL ELEQ,TEQ,TEQM,TENDD,DTD,DXD,DYD,TDECAY,WSCALE
      REAL RFRIC,STRESS0
      real A,B,C,phik,y
      real q0,qm,X,F,Fk,Fm,Cxm,Cym
      real XX,YY,XXT,YYT,UB,V,HB,BigX
      real ysouth,ynorth,xwest
c     
c Kelvin Wave Forcing 
c
      cs=0.5*dt/dx

      Fk(2)=cs*q0(1)+q0(2)+dt*Cxm(0,2)

      do i=3,nx
         
         Fk(i)=q0(i)+dt*Cxm(0,i)

      enddo
c
c Rossby Waves Forcing
c
      do m=1,modes-1

         Fm(m+1,1)=qm(m+1,1)+dt/(2.0*m+1.0)*(
     &        m*Cxm(m+1,1)-sqrt(m*(m+1.0))*Cxm(m-1,1)+
     &        sqrt(2.0*(m+1.0))*(rfric*Cym(m,1)-
     &        (Cym(m,2)-Cym(m,1))/dx))
         
         do i=2,nx-2
            
            Fm(m+1,i)=qm(m+1,i)+dt/(2.0*m+1.0)*(
     &           m*Cxm(m+1,i)-sqrt(m*(m+1.0))*Cxm(m-1,i)+
     &           sqrt(2.0*(m+1.0))*(rfric*Cym(m,i)-
     &           (Cym(m,i+1)-Cym(m,i-1))/(2.0*dx)))

         enddo

         Fm(m+1,nx-1)=cs*qm(m+1,nx)/(2.0*m+1.0)+
     &        qm(m+1,nx-1)+dt/(2.0*m+1.0)*(
     &        m*Cxm(m+1,nx-1)-sqrt(m*(m+1.0))*Cxm(m-1,nx-1)+
     &        sqrt(2.0*(m+1.0))*(rfric*Cym(m,nx-1)-
     &        (Cym(m,nx)-Cym(m,nx-2))/(2.0*dx)))

      enddo
c
c Putting all the equations in one vector (Vector F)
c
      call VectorF

      return
      end

c*********************************************************
c Subroutine to generate the vector X. Vector X contains *
c the kelvin wave (q0) and all the rossby waves (qm).    * 
c*********************************************************

      subroutine VectorX

      include 'zeq.common'

c      implicit none
      integer i,j,m
      INTEGER NX,NY,NXP,NYP,IALPHA,ISTART,NT,NTIMES,NPRINT
      INTEGER NSEG,NTAPE,NREWND,NSTART,NATM,NATMR,MODES,NDIM
      REAL DT,DX,DY,TPLMIN,HK,PI,SQRT2,TD,TZERO,XE,YN,YS
      REAL OMEGA,ENDT,RALPHA,HEQUIV,XWD,XED,YSD,YND,CWAVE
      REAL ELEQ,TEQ,TEQM,TENDD,DTD,DXD,DYD,TDECAY,WSCALE
      REAL RFRIC,STRESS0
      real A,B,C,phik,y
      real q0,qm,X,F,Fk,Fm,Cxm,Cym
      real XX,YY,XXT,YYT,UB,V,HB,BigX
      real ysouth,ynorth,xwest
c
c Components of X related to the Kelvin wave 
c
      do i=1,nx-1
         
         X(i)=q0(i+1)

      enddo
c
c Components of X related to the Rossby waves 
c
      do m=1,modes-1
         do i=1,nx-1
            
            j=(nx-1)*m+i
            X(j)=qm(m+1,i)

         enddo
      enddo

      return
      end

c*********************************************************
c Subroutine to generate the vector F. Vector F contains *
c the tendencies of vector X.  In other words, F contains*
c the RHS of the equations.                              *
c*********************************************************

      subroutine VectorF

      include 'zeq.common'

c      implicit none
      integer i,j,m
      INTEGER NX,NY,NXP,NYP,IALPHA,ISTART,NT,NTIMES,NPRINT
      INTEGER NSEG,NTAPE,NREWND,NSTART,NATM,NATMR,MODES,NDIM
      REAL DT,DX,DY,TPLMIN,HK,PI,SQRT2,TD,TZERO,XE,YN,YS
      REAL OMEGA,ENDT,RALPHA,HEQUIV,XWD,XED,YSD,YND,CWAVE
      REAL ELEQ,TEQ,TEQM,TENDD,DTD,DXD,DYD,TDECAY,WSCALE
      REAL RFRIC,STRESS0
      real A,B,C,phik,y
      real q0,qm,X,F,Fk,Fm,Cxm,Cym
      real XX,YY,XXT,YYT,UB,V,HB,BigX
      real ysouth,ynorth,xwest
c
c Forcing for the Kelvin wave 
c
      do i=1,nx-1

         F(i)=Fk(i+1)

      enddo
c
c Forcing for the Rossby waves 
c
      do m=1,modes-1
         do i=1,nx-1
            
            j=(nx-1)*m+i
            F(j)=Fm(m+1,i)

         enddo
      enddo

      return
      end

c**************************************
c Subroutine to transforms the vector * 
c X back into the fields q0 and qm.   *
c**************************************

      subroutine VectorQ

      include 'zeq.common'

c      implicit none
      integer i,j,m
      INTEGER NX,NY,NXP,NYP,IALPHA,ISTART,NT,NTIMES,NPRINT
      INTEGER NSEG,NTAPE,NREWND,NSTART,NATM,NATMR,MODES,NDIM
      REAL DT,DX,DY,TPLMIN,HK,PI,SQRT2,TD,TZERO,XE,YN,YS
      REAL OMEGA,ENDT,RALPHA,HEQUIV,XWD,XED,YSD,YND,CWAVE
      REAL ELEQ,TEQ,TEQM,TENDD,DTD,DXD,DYD,TDECAY,WSCALE
      REAL RFRIC,STRESS0
      real A,B,C,phik,y
      real q0,qm,X,F,Fk,Fm,Cxm,Cym
      real XX,YY,XXT,YYT,UB,V,HB,BigX
      real ysouth,ynorth,xwest
c
c The Kelvin wave 
c
      do i=1,nx-1
         
         q0(i+1)=X(i)

      enddo
c
c The Rossby waves 
c
      do i=1,nx-1
         
         qm(1,i)=sqrt(2.0)*Cym(0,i)

      enddo

      do m=1,modes-1
         do i=1,nx-1
            
            j=(nx-1)*m+i
            qm(m+1,i)=X(j)

         enddo
      enddo
c
c Western Boundary Condition
c
      q0(1)=0.707*qm(2,1)+0.176*qm(4,1)+0.095*qm(6,1)+
     &     0.063*qm(8,1)+0.046*qm(10,1)+0.036*qm(12,1)+
     &     0.029*qm(14,1)+0.024*qm(16,1)+0.021*qm(18,1)+
     &     0.018*qm(20,1)+0.016*qm(22,1)+0.014*qm(24,1)+   
     &     0.013*qm(26,1)+0.012*qm(28,1)+0.011*qm(30,1) 
c     
c Eastern Boundary Condition
c
      qm(1,nx)=sqrt(2.0)*Cym(0,nx)
      qm(2,nx)=sqrt(0.5)*q0(nx)+Cym(1,nx)
      
      do m=2,modes-1
      
         qm(m+1,nx)=sqrt(m/(m+1.0))*qm(m-1,nx)+
     &        sqrt(2.0/(m+1.0))*Cym(m,nx)
      
      enddo

      return
      end

***********************************************
* FINITE DIFFERENCE PDE INTEGRATION SUBROUTINE*
***********************************************

      subroutine Matrix

      include 'zeq.common'

c      implicit none
      integer i,j,k,m
      real cs
      INTEGER NX,NY,NXP,NYP,IALPHA,ISTART,NT,NTIMES,NPRINT
      INTEGER NSEG,NTAPE,NREWND,NSTART,NATM,NATMR,MODES,NDIM
      REAL DT,DX,DY,TPLMIN,HK,PI,SQRT2,TD,TZERO,XE,YN,YS
      REAL OMEGA,ENDT,RALPHA,HEQUIV,XWD,XED,YSD,YND,CWAVE
      REAL ELEQ,TEQ,TEQM,TENDD,DTD,DXD,DYD,TDECAY,WSCALE
      REAL RFRIC,STRESS0
      real A,B,C,phik,y
      real q0,qm,X,F,Fk,Fm,Cxm,Cym
      real XX,YY,XXT,YYT,UB,V,HB,BigX
      real ysouth,ynorth,xwest

      cs=0.5*dt/dx
c
c Lower Diagonal (matrix A)
c
      do i=2,nx-2
         
         A(i)=-cs    
         
      enddo

      A(nx-1)=-2.0*cs

      do m=1, modes-1

         k=(nx-1)*m+1
         A(k)=0.0

         do i=2,nx-1

            j=(nx-1)*m+i
            A(j)=cs/(2.0*m+1.0)

         enddo

      enddo
c
c Main Diagonal (matrix B)
c
      do i=1,nx-2
        
         B(i)=(1.0+rfric*dt)    
         
      enddo

      B(nx-1)=(1.0+rfric*dt+2.0*cs)

      do m=1, modes-1

         j=(nx-1)*m+1
         B(j)=1.0+rfric*dt+2.0*cs/(2.0*m+1.0)

         do i=1,nx-2
            
            k=j+i
            B(k)=(1.0+rfric*dt)

         enddo

         B(k+1)=1.0+rfric*dt+2.0*cs/(2.0*m+1.0)

      enddo
c
c Upper Diagonal
c      
      do i=1,nx-2
         
         C(i)=cs            

      enddo

      do m=1, modes-1

         k=(nx-1)*m
         C(k)=0.0
         C(k+1)=-2.0*cs/(2.0*m+1.0)

         do i=2,nx-2
            
            j=(nx-1)*m+i
            C(j)=-cs/(2.0*m+1.0)
            
         enddo
      enddo

      return
      end

*********************************
* TRIDIAGONAL SOLVING SUBROUTINE*
*********************************
c
c Check numerical recipes 
c
      subroutine tridiagonal

      include 'zeq.common'

c      implicit none
      integer n,j 
      real gam(ndim),bet
      INTEGER NX,NY,NXP,NYP,IALPHA,ISTART,NT,NTIMES,NPRINT
      INTEGER NSEG,NTAPE,NREWND,NSTART,NATM,NATMR,MODES,NDIM
      REAL DT,DX,DY,TPLMIN,HK,PI,SQRT2,TD,TZERO,XE,YN,YS
      REAL OMEGA,ENDT,RALPHA,HEQUIV,XWD,XED,YSD,YND,CWAVE
      REAL ELEQ,TEQ,TEQM,TENDD,DTD,DXD,DYD,TDECAY,WSCALE
      REAL RFRIC,STRESS0
      real A,B,C,phik,y
      real q0,qm,X,F,Fk,Fm,Cxm,Cym
      real XX,YY,XXT,YYT,UB,V,HB,BigX
      real ysouth,ynorth,xwest
c
c Reduce matrix, which is store by bands
c
      n=ndim
      If (B(1) .eq. 0.0) pause
      bet=B(1)
      X(1)=F(1)/bet

      do j=2,n
         gam(j)=C(j-1)/bet
         bet=B(j)-A(j)*gam(j)
         If (bet .eq. 0.0) pause
         X(j)=(F(j)-A(j)*X(j-1))/bet
      enddo
c
c Back substitute. 
c
      do j=n-1,1,-1
         X(j)=X(j)-gam(j+1)*X(j+1)
      enddo

      return
      end


EOF

rm -f axcforce.f
cat>axcforce.f<<EOF
       SUBROUTINE CFORCE
 
c this routine supplies wind stress forcing to the ocean dynamics model
c For coupled ocean/atmosphere model the wind stress is derived from
c the atmosphere model, and is saved in the array HTAU.

      INCLUDE 'zeq.common'
 
      DIMENSION  HT(34,30),XOFF(2),YOFF(2)

      COMMON/ZDATA/WM(30,34,2,12),DIVM(30,34,12),Q0O(30,34),
     A SSTM(30,34,12),UO(30,34),VO(30,34),DO(30,34),HTAU(34,30,2),
     B WEM(30,34,12),TO(30,34),UV(30,34,2,12),US(30,34),VS(30,34),
     C WP(30,34),DT1(30,34,7),TT(30,34),UV1(30,34),UV2(30,34),
     D WM1(30,34),UAT(30,34),VAT(30,34),DIVT(30,34),HZC(30,34,12),
     F HLF(30,34,12),TAUMX(1020,1020),TAUMY(1020,1020)

c following are the wind field parameters
      DATA   WXW/+101.25/,WXE/+286.875/,WYS/-29./,WYN/+29./
      DATA   IW/34/,JW/30/

c  offsets for the staggering of F and G
      DATA XOFF/0.5E0,1.E0/,YOFF/0.5E0,1.E0/
c......................................................................
c some initialization computations. Gridding factors...

      WDX=(WXE-WXW)/(IW-1)
      AWX=1.E0+(XWD-WXW)/WDX
      WDXS=ELEQ*DX/WDX
      WDY=(WYN-WYS)/(JW-1)
      AWY=1.E0+(YSD-WYS)/WDY
      WDYS=ELEQ*DY/WDY
 
c on the initial time step, only G is needed, unless restarting.
 
      ISM=1
      IEM=NX
      JSM=1
      JEM=NY
      L1=1
      L2=2
 
c get the wind stress by accessing the atmosphere model array HTAU...
c the stress data are interpolated onto ocean grids...
c the code assumes the wind data domain includes the entire ocean domain

      DO 130 L=L1,L2
      ASSIGN 119 TO LFG
      IF(L.EQ.2) ASSIGN 118 TO LFG
      AAWX=AWX-XOFF(L)*WDXS
      AAWY=AWY-YOFF(L)*WDYS
 
      DO 110 J=1,JW
      J13=JW+1-J
      DO 111 I=1, IW
      HT(I,J)=HTAU(I,J,L)
111   CONTINUE
110   CONTINUE
 
      IEG = IEM + (L-1)

      DO 125 I=ISM,IEG
      RX=AAWX+WDXS*I
      IX=RX
      RX=RX-IX

c ocean basin segment boundary check

      JEG = JEM+L-1

      DO 120 J=1,JEG
      RY=AAWY+WDYS*J
      JY=RY
      RY=RY-JY
      FA=HT(IX,JY)+RY*(HT(IX ,JY+1)-HT(IX ,JY))
      FB=HT(IX+1,JY)+RY*(HT(IX+1,JY+1)-HT(IX+1,JY))
      FA = FA + RX * ( FB - FA )
      GO TO LFG,(118,119)
118   YYT(J,I)=FA
      GO TO 120
119   XXT(J,I)=FA
      
120   CONTINUE

125   CONTINUE
130   CONTINUE

      do i=1,nx
         do j=1,ny
               
            XX(i,j)=wscale*XXT(j,i)
            YY(i,j)=wscale*YYT(j,i)
               
         enddo
      enddo

      RETURN
      END

EOF

rm -f axatmos.f
cat>axatmos.f<<EOF
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


EOF

rm -f axdridag.f
cat>axtridag.f<<EOF
      SUBROUTINE TRIDAG(A,B,E,D,X,JS,N)
 
c a tridiagonal linear system solver for the system:
c   A(I)*X(I-1)+B(I)*X(I)+C(I)*X(I+1)=D(I)...I=JS,N  A(JS)=0,C(N)=0
c   ARRAY D IS DESTROYED IN THE COMPUTATION
c   ASSUMES C(J) = A(J+1)
c   E IS TEMPORARY STORAGE
 
      DIMENSION A(N+1),B(N+1),E(N+1),D(N+1),X(N+1)
 
      E(JS+1) = A(JS+2)/B(JS+1)
      D(JS+1)=D(JS+1)/B(JS+1)
 
      JN = JS+2
      DO 10 I=JN,N
        DN=B(I)-A(I)*E(I-1)
        E(I)=A(I+1)/DN
        D(I)=(D(I)-A(I)*D(I-1))/DN
10    CONTINUE
 
      X(N)=D(N)
      I=N
      DO 20 II=JN,N
         I=I-1
         X(I)=D(I)-E(I)*X(I+1)
20    CONTINUE
 
      RETURN
      END

EOF

rm -f axinitdat_FOM.f
cat>axinitdat_FOM.f<<EOF
       SUBROUTINE INITDAT
 
c this routine reads in and sets initial values of data arrays for
c the coupled model.

      INCLUDE 'zeq.common'
 
      COMMON/ZDATA/WM(30,34,2,12),DIVM(30,34,12),Q0O(30,34),
     A SSTM(30,34,12),UO(30,34),VO(30,34),DO(30,34),HTAU(34,30,2),
     B WEM(30,34,12),TO(30,34),UV(30,34,2,12),US(30,34),VS(30,34),
     C WP(30,34),DT1(30,34,7),TT(30,34),UV1(30,34),UV2(30,34),
     D WM1(30,34),UAT(30,34),VAT(30,34),DIVT(30,34),HZC(30,34,12),
     F HLF(30,34,12),TAUMX(1020,1020),TAUMY(1020,1020)

c read in the mean wind, sst, divergence, and currents data
 
      DIMENSION XX1(30,34),XX2(30,34),XX3(30,34,2),
     &       XX4(30,34),XX5(30,34,2),XX6(30,34)

      DIMENSION XX7(1020,1020),XX8(1020,1020)

      DIMENSION YY1(30,34,12),YY2(30,34,12),YY4(30,34,12)
      DIMENSION YY3(30,34,2,12),YY5(30,34,2,12)

      READ(53,555) ((XX7(i,j),j=1,1020),i=1020,1,-1)
      READ(54,555) ((XX8(i,j),j=1,1020),i=1020,1,-1)
555   FORMAT(6E12.4)

      do k2=1,1020
         do k1=1,1020

            TAUMX(k1,k2)=XX7(k1,k2)
            TAUMY(k1,k2)=XX8(k1,k2)
               
         enddo
      enddo

c
C********************************************************
C Annual Mean Basic State
C
c      READ(84,REC=1) ((XX1(I,J),J=1,34),I=30,1,-1)           ! From BSTATES
c      READ(87,REC=1) ((XX4(I,J),J=1,34),I=30,1,-1)           ! From BSTATES
c      READ(89,REC=1) (((XX5(I,J,K),J=1,34),I=30,1,-1),K=1,2) ! From BSTATES
c      READ(80,REC=1) ((XX6(I,J),J=1,34),I=30,1,-1)            ! From BSTATES
c
      DO 303 L=1,12
c
c      READ(83) ((XX2(I,J),I=1, 30),J=1,34)
c      READ(82) (((XX3(I,J,K),I=1, 30),J=1, 34),K=1,2)
c
c      READ(80,REC=1) ((XX6(I,J),J=1,34),I=30,1,-1)
      READ(84,555) ((XX1(I,J),I=1, 30),J=1, 34)
          

      write(6,*) 'stop'
*      stop

      READ(83,555) ((XX2(I,J),I=1, 30),J=1, 34)
      READ(82,555) (((XX3(I,J,K),I=1, 30),J=1, 34),K=1,2)
      READ(87,555) ((XX4(I,J),I=1, 30), J=1, 34)
      READ(89,555) (((XX5(I,J,K),I=1, 30),J=1, 34),K=1,2)
c
      DO I=1,30
      DO J=1,34
      SSTM(I,J,L)=XX1(I,J)
      DIVM(I,J,L)=XX2(I,J)
      WEM(I,J,L)=XX4(I,J)
c      HLF(I,J,L)=XX6(I,J)         ! h at each mean basic state
      DO K=1,2
      WM(I,J,K,L)=XX3(I,J,K)
      UV(I,J,K,L)=XX5(I,J,K)
      ENDDO
      ENDDO
      ENDDO
303   CONTINUE
c
C********************************************************
C Annual Mean Basic State
c
      DO 307 L=1,12
      DO I=1,30
      DO J=1,34

      YY1(I,J,L)=(SSTM(I,J,1)+SSTM(I,J,2)+SSTM(I,J,3)+
     &        SSTM(I,J,4)+SSTM(I,J,5)+SSTM(I,J,6)+
     &        SSTM(I,J,7)+SSTM(I,J,8)+SSTM(I,J,9)+
     &        SSTM(I,J,10)+SSTM(I,J,11)+SSTM(I,J,12))/12.0
      
      YY2(I,J,L)=(DIVM(I,J,1)+DIVM(I,J,2)+DIVM(I,J,3)+
     &        DIVM(I,J,4)+DIVM(I,J,5)+DIVM(I,J,6)+
     &        DIVM(I,J,7)+DIVM(I,J,8)+DIVM(I,J,9)+
     &        DIVM(I,J,10)+DIVM(I,J,11)+DIVM(I,J,12))/12.0

      YY4(I,J,L)=(WEM(I,J,1)+WEM(I,J,2)+WEM(I,J,3)+
     &        WEM(I,J,4)+WEM(I,J,5)+WEM(I,J,6)+
     &        WEM(I,J,7)+WEM(I,J,8)+WEM(I,J,9)+
     &        WEM(I,J,10)+WEM(I,J,11)+WEM(I,J,12))/12.0

      DO K=1,2

      YY3(I,J,K,L)=(WM(i,J,K,1)+WM(i,J,K,2)+WM(i,J,K,3)+
     &        WM(i,J,K,4)+WM(i,J,K,5)+WM(i,J,K,6)+
     &        WM(i,J,K,7)+WM(i,J,K,8)+WM(i,J,K,9)+
     &        WM(i,J,K,10)+WM(i,J,K,11)+WM(i,J,K,12))/12.0

      YY5(I,J,K,L)=(UV(I,J,K,1)+UV(I,J,K,2)+UV(I,J,K,3)+
     &        UV(I,J,K,4)+UV(I,J,K,5)+UV(I,J,K,6)+
     &        UV(I,J,K,7)+UV(I,J,K,8)+UV(I,J,K,9)+
     &        UV(I,J,K,10)+UV(I,J,K,11)+UV(I,J,K,12))/12.0

      ENDDO
      ENDDO
      ENDDO
 307  CONTINUE






      DO L=1,12
         DO I=1,30
            DO J=1,34

               SSTM(I,J,L)=YY1(I,J,L)
               DIVM(I,J,L)=YY2(I,J,L)
               WEM(I,J,L)=YY4(I,J,L)
               
               DO K=1,2

                  WM(I,J,K,L)=sqrt(1.0)*YY3(I,J,K,L)
                  UV(I,J,K,L)=YY5(I,J,K,L)

               ENDDO
            ENDDO
         ENDDO
      ENDDO
c
c***********************************************************
c
c if not restarting, zero out arrays q0o,uo,vo,do,to
c
C      IF (NSTART .GT. 0 .AND. NATMR .GT. 0) GO TO 311
 
      DO 310 I=1, 30
      DO 310 J=1, 34
      Q0O(I,J)=0.E0
      UO(I,J)=0.E0
      VO(I,J)=0.E0
      DO(I,J)=0.E0
      TO(I,J)=0.E0
310   CONTINUE
311   CONTINUE

c
c The perturbation for T0 has to be done only in the atmos grid,
c otherwise the lon-lat plots of T0 presents the frame at the 
c boundaries.
c
      DO 350 I=6,25
      DO 351 J=6,32
C      TO(I,J)=0.0
 351  CONTINUE
 350  CONTINUE

      do i=1,1439
C         BigX(i)=0.0
      enddo

      do i=1,899
C         X(i)=BigX(i)
      enddo

      RETURN
      END


EOF

rm -f axhermite.f
cat>axhermite.f<<EOF
c*******************************
c Subroutine Hermite functions *
c*******************************

      subroutine hermite

      include 'zeq.common'

      integer m,j
      real sum1(modes-1),sum2(modes-1)
      real enorm

      enorm=sqrt(pi)
c
c Hermite functions for m=0 and m=1
c
      do j=1,ny
         
         phik(0,j)=sqrt(1.0/enorm)*exp(-0.5*y(j)*y(j))
         phik(1,j)=sqrt(2.0/enorm)*y(j)*
     &        exp(-0.5*y(j)*y(j))

      enddo
c
c Hermite functions for m=2,3....m=modes
c         
      do m=1,modes-1
         do j=1,ny
            
            phik(m+1,j)=sqrt(2.0/(m+1.0))*y(j)*phik(m,j)-
     &           sqrt(m/(m+1.0))*phik(m-1,j)

         enddo
      enddo
c
c Checking normalization of the Hermite functions
c
      do m=1,modes-1

         sum1(m)=0.0
         sum2(m)=0.0

         do j=1,ny

            sum1(m)=sum1(m)+phik(m,j)*phik(m,j)*dy
            sum2(m)=sum2(m)+phik(m,j)*phik(m+1,j)*dy            

         enddo
      enddo
      
      return
      end

EOF

rm -f xmym.f
cat>axxmym.f<<EOF
c*******************************************
c Subroutine to transform the stress field *
c into the coefficients Xm and Ym          *
c*******************************************

      subroutine XmYm

      include 'zeq.common'

      integer i,j,m
c
c Computing the stress coefficients Xm and Ym
c
      do m=1,modes+1      
         do i=1,nx

            Cxm(m-1,i)=0.0
            Cym(m-1,i)=0.0

            do j=1,ny
               
               Cxm(m-1,i)=Cxm(m-1,i)+XX(i,j)*phik(m-1,j)*dy               
               Cym(m-1,i)=Cym(m-1,i)+YY(i,j)*phik(m-1,j)*dy

            enddo
         enddo
      enddo

      return
      end

EOF

rm -f axcurrents.f
cat>axcurrents.f<<EOF
c*********************
c Compute the fields *
c*********************

      subroutine currents

      include 'zeq.common'

      integer i,j,m
      real rm(0:modes-2,nx),vm(0:modes-1,nx-1)
      real qxyt(0:modes,nx,ny)
      real Q(nx,ny),P(nx,ny),R(nx,ny)
      real U(nx,ny),VV(nx-1,ny),h(nx,ny)
c
c Computing the total Q and R.  U and h have different
c truncation, that explains the definition of Q and P.
c
      do i=1,nx
         do j=1,ny
               
            qxyt(0,i,j)=phik(0,j)*q0(i)               

         enddo
      enddo

      do m=1,modes
         do i=1,nx
            do j=1,ny
      
               qxyt(m,i,j)=phik(m,j)*qm(m,i)

            enddo
         enddo
      enddo

      do i=1,nx
         do j=1,ny
            Q(i,j)=0.0
            do m=0,modes

               Q(i,j)=Q(i,j)+qxyt(m,i,j)
               
            enddo
         enddo
      enddo

      do i=1,nx
         do j=1,ny
            P(i,j)=0.0
            do m=0,modes-2

               P(i,j)=P(i,j)+qxyt(m,i,j)
               
            enddo
         enddo
      enddo

c
c Compute the diagnostic rm and the total R. Use
c the eastern boundary condition for rm. 
c
      do m=1,modes-1
         do i=1,nx
            
            rm(m-1,i)=sqrt((m+1.0)/m)*qm(m+1,i)
     &           -sqrt(2.0/m)*Cym(m,i)
            
         enddo
      enddo

      do i=1,nx
         do j=1,ny
            R(i,j)=0.0
            do m=0,modes-2
         
               R(i,j)=R(i,j)+phik(m,j)*rm(m,i)
               
            enddo
         enddo
      enddo
c
c Computing the fields u and h
c
      do i=1,nx
         do j=1,ny
            
            U(i,j)=0.5*(P(i,j)-R(i,j))
            h(i,j)=0.5*(Q(i,j)+R(i,j))

         enddo
      enddo
c
c Computing the diagnostic vm and the total V
c
      vm(0,1)=sqrt(0.5)*(2.0*(qm(1,2)-qm(1,1))/dx-
     &     Cxm(1,1)+sqrt(2.0)*(rfric*Cym(0,1)-
     &     (Cym(0,2)-Cym(0,1))/dx))

      do i=2,nx-1
         
         vm(0,i)=sqrt(0.5)*((qm(1,i+1)-qm(1,i-1))/dx-
     &        Cxm(1,i)+sqrt(2.0)*(rfric*Cym(0,i)-
     &        (Cym(0,i+1)-Cym(0,i-1))/(2.0*dx)))

      enddo

      do m=1,modes-1

         vm(m,1)=sqrt((m+1.0)/2.0)/(2.0*m+1.0)*
     &        (2.0*(qm(m+1,2)-qm(m+1,1))/dx-
     &        Cxm(m+1,1)-sqrt(m/(m+1.0))*Cxm(m-1,1)+
     &        sqrt(2.0/(m+1.0))*(rfric*Cym(m,1)-
     &        (Cym(m,2)-Cym(m,1))/dx))

         do i=2,nx-1
            
            vm(m,i)=sqrt((m+1.0)/2.0)/(2.0*m+1.0)*
     &           ((qm(m+1,i+1)-qm(m+1,i-1))/dx-
     &           Cxm(m+1,i)-sqrt(m/(m+1.0))*Cxm(m-1,i)+
     &           sqrt(2.0/(m+1.0))*(rfric*Cym(m,i)-
     &           (Cym(m,i+1)-Cym(m,i-1))/(2.0*dx)))

         enddo
      enddo

      do i=1,nx-1
         do j=1,ny
            VV(i,j)=0.0
            do m=0,modes-1

               VV(i,j)=VV(i,j)+phik(m,j)*vm(m,i)

            enddo
         enddo
      enddo
c
c Transforming the arrays U, VV and h into UB,V and HB.
c These arrays are the ones to be used in the sst 
c subroutine in ZC code.  The coupling requires HB in mts
c (so we multiplying h by 150.0) and UB ad V must be
c non-dimensional. Since, the currents are leaving this
c code non-dimensional, they have to be multiply by 2.5
c in the main program to be plotted. In other words the'
c DUMMY array has to be multiply by 2.5m/sec to get the
c currents in m/sec. 
c
      do i=1,nx
         do j=1,ny 

            UB(j,i)=U(i,j)
            HB(j,i)=150.0*h(i,j)

         enddo
      enddo

      do i=1,nx-1
         do j=1,ny 

            V(j,i)=VV(i,j)

         enddo
      enddo

      return
      end

EOF

rm -f axsetup2.f
cat>axsetup2.f<<EOF
      SUBROUTINE SETUP2

c routine for reading or initializing variables related to SST physics
c and atmospheric physics....

      COMMON/SGRID/XWRS,DXRS,YNRS,DYRS,NYPS

      COMMON/SSTPAR/GAM1,GAM2,TDA1,TDB1,TDA2,TDB2,TLOSS

      COMMON/ATMPAR/ALPHA,EPS,BETA,CTOL,IMIN,IMAX,ISTEP

      COMMON/INIT/IC2

      COMMON/STRPAR/TCUT,TSEED,TAMP1,TAMP2,JLEFT,JRT,ITOP,IBOT,NIC

      LOGICAL ISTEP

      SAVE

c  define ssta and atmosphere grid for plotting package

      XWRS=101.25
      DXRS=5.625
      YNRS=29.0
      DYRS=2.0
      NYPS=30

c read in init parameter
      READ (5,2) IC2
 
c read in ssta parameters......
      READ (5,3) GAM1,GAM2,TDA1,TDB1,TDA2,TDB2,TLOSS

c read in atmosphere model parameters
      READ (5,4) ALPHA,EPS,BETA,CTOL,IMIN,IMAX,ISTEP

c read in external wind stress parameters (for initializing coupled model
c or adding intraseasonal noise..

      READ (5,5) TCUT,TSEED,TAMP1,TAMP2,JLEFT,JRT,ITOP,IBOT,NIC

2     FORMAT(8X,I15)
3     FORMAT(8X,F15.0)
4     FORMAT(4(8X,F15.0 / ),2(8X,I15 / ),8X,L15 )
5     FORMAT(4(8X,F15.0 / ),4(8X,I15 / ),8X,I15 )
   

      WRITE(6,6) IC2 
      WRITE(6,7) GAM1,GAM2,TDA1,TDB1,TDA2,TDB2,TLOSS
      WRITE(6,8) ALPHA,EPS,BETA,CTOL,IMIN,IMAX,ISTEP
      WRITE(6,9) TCUT,TSEED,TAMP1,TAMP2,JLEFT,JRT,ITOP,IBOT,NIC

6     FORMAT ( ' IC2    =',I15   //)

7     FORMAT ( ' GAM1   =',E15.7, 5X, ' GAM2   =',E15.7 /
     3         ' TDA1   =',E15.7, 5X, ' TDB1   =',E15.7 /
     5         ' TDA2   =',E15.7, 5X, ' TDB2   =',E15.7 /
     7         ' TLOSS  =',E15.7 // )

8     FORMAT ( ' ALPHA  =',E15.7, 5X, ' EPS    =',E15.7 /
     A         ' BETA   =',E15.7, 5X, ' CTOL   =',E15.7 /
     B         ' IMIN   =',I15  , 5X, ' IMAX   =',I15   /
     B         ' ISTEP  =',L15   // )

9     FORMAT ( ' TCUT   =',E15.7, 5X, ' TSEED  =',E15.7 /
     F         ' TAMP1  =',E15.7, 5X, ' TAMP2  =',E15.7 /
     F         ' JLEFT  =',I15,   5X, ' JRT    =',I15   /
     F         ' ITOP   =',I15,   5X, ' IBOT   =',I15   /
     F         ' NIC    =',I15   // )

      RETURN
      END

EOF

rm -f axsst.f
cat>axsst.f<<EOF
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


EOF

rm -f axfft2c.f
cat>axfft2c.f<<EOF
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

EOF

 
foreach file (main_cc openfl_FOM taumatrix  setup constc waves-1 cforce atmos tridag initdat_FOM  hermite xmym currents setup2  sst fft2c)

${FC} -c ax${file}.f

end
rm -f zcout
${FC} -o zcout axmain_cc.o axopenfl_FOM.o axtaumatrix.o axsetup.o axconstc.o axwaves-1.o axcforce.o axatmos.o \
       axtridag.o axinitdat_FOM.o axhermite.o axxmym.o axcurrents.o axsetup2.o  axsst.o axfft2c.o 

./zcout


cp nino.gdat nino_mu${MU}_gam${GAM}.gdat


rm -f plo.gnu 
cat>plo.m<<EOF
load nino.gdat
a=size(nino);
LL=[1:a]/12
plot(LL,nino,'Linewidth',1)
grid
xlabel('time (years)')
ylabel('Nino 3 SSTA')
print -depsc nino3.ps
EOF
matlab<<EOF
plo
EOF
open nino3.ps

end

exit
