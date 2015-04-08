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

