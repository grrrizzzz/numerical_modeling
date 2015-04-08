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

