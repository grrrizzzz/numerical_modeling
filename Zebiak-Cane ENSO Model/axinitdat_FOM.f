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


