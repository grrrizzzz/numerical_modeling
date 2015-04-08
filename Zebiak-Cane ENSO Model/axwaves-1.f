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


