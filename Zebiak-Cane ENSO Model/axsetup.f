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

