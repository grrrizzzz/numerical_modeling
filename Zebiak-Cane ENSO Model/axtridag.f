      SUBROUTINE TRIDAG(A,B,E,D,X,JS,N)
 
c a tridiagonal linear system solver for the system:
c   A(I)*X(I-1)+B(I)*X(I)+C(I)*X(I+1)=D(I)...I=JS,N  A(JS)=0,C(N)=0
c   ARRAY D IS DESTROYED IN THE COMPUTATION
c   ASSUMES C(J) = A(J+1)
c   E IS TEMPORARY STORAGE
 
      DIMENSION A(N+1),B(N+1),E(N+1),D(N+1),X(N+1)
 
      E(JS+1) = A(JS+2)/B(JS+1)
      D(JS+1)=D(JS+1)/B(JS+1)
 
      JN = JS+2
      DO 10 I=JN,N
        DN=B(I)-A(I)*E(I-1)
        E(I)=A(I+1)/DN
        D(I)=(D(I)-A(I)*D(I-1))/DN
10    CONTINUE
 
      X(N)=D(N)
      I=N
      DO 20 II=JN,N
         I=I-1
         X(I)=D(I)-E(I)*X(I+1)
20    CONTINUE
 
      RETURN
      END

