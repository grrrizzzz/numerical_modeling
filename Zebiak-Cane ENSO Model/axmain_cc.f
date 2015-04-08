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
      FATM = 0.95
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
 

