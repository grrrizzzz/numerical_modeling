c*******************************
c Subroutine Hermite functions *
c*******************************

      subroutine hermite

      include 'zeq.common'

      integer m,j
      real sum1(modes-1),sum2(modes-1)
      real enorm

      enorm=sqrt(pi)
c
c Hermite functions for m=0 and m=1
c
      do j=1,ny
         
         phik(0,j)=sqrt(1.0/enorm)*exp(-0.5*y(j)*y(j))
         phik(1,j)=sqrt(2.0/enorm)*y(j)*
     &        exp(-0.5*y(j)*y(j))

      enddo
c
c Hermite functions for m=2,3....m=modes
c         
      do m=1,modes-1
         do j=1,ny
            
            phik(m+1,j)=sqrt(2.0/(m+1.0))*y(j)*phik(m,j)-
     &           sqrt(m/(m+1.0))*phik(m-1,j)

         enddo
      enddo
c
c Checking normalization of the Hermite functions
c
      do m=1,modes-1

         sum1(m)=0.0
         sum2(m)=0.0

         do j=1,ny

            sum1(m)=sum1(m)+phik(m,j)*phik(m,j)*dy
            sum2(m)=sum2(m)+phik(m,j)*phik(m+1,j)*dy            

         enddo
      enddo
      
      return
      end

