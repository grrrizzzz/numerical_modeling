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

