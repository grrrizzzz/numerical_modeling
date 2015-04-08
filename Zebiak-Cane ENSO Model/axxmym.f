c*******************************************
c Subroutine to transform the stress field *
c into the coefficients Xm and Ym          *
c*******************************************

      subroutine XmYm

      include 'zeq.common'

      integer i,j,m
c
c Computing the stress coefficients Xm and Ym
c
      do m=1,modes+1      
         do i=1,nx

            Cxm(m-1,i)=0.0
            Cym(m-1,i)=0.0

            do j=1,ny
               
               Cxm(m-1,i)=Cxm(m-1,i)+XX(i,j)*phik(m-1,j)*dy               
               Cym(m-1,i)=Cym(m-1,i)+YY(i,j)*phik(m-1,j)*dy

            enddo
         enddo
      enddo

      return
      end

