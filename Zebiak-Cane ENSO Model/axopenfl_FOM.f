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

