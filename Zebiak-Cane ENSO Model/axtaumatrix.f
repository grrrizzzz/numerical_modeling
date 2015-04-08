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
 

