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

