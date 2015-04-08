close
clear
clf

N=300;
M=300;
%TI=20000;
TI=1000;
DX=25*1E3;
DY=DX;
DT=60;
DPHI=DY/40e6*360;

RHO=1000;


LATS=10;
LON=DX/40e6*360*[1:N];

LAT=LATS+DPHI*[1:M];
YS=LATS*40000e3/360;
Y=YS+DX*[1:M];
F=2*2*pi/(86400)*sin(LAT/180*pi);
BETA=F./Y;

OMEGA = 1.6;


%%%%%%%%%% INITIALIZE

H=ones(M,N)*4e3;
PSI=ones(M,N);
PSI(:,end)=LON';
PSI(:,1)=LON';
PSI(end,:)=LAT;
PSI(1,:)=LAT;
ZET=zeros(M,N);
ZET(:,end)=LON';
ZET(:,1)=LON';
ZET(end,:)=LAT;
ZET(1,:)=LAT;
TAUY=zeros(M,N);
NONL1=zeros(M,N);
NONL=zeros(M,N);
CURL=zeros(M,N);
DISS=zeros(M,N);
DIV1=zeros(M,N);
DIV=zeros(M,N);
ZN=zeros(M,N);
PHIN=zeros(M,N);
PSINEW=zeros(M,N);
ZETNEW=zeros(M,N);
PSIOLD=zeros(M,N);
ZETOLD=zeros(M,N);


for i=1:M
    for j=1:N
TAUX(i,j)=0.01.*cos(2*pi*LAT(j)/360/0.2);
    end
end


AH=3e-8;

eps_sor=0.001;
anormf=0;
anorm=1;
itt=0;

for itime=1:TI
    
    if (itime==1)
        
        for I=2:M-1
        for J=2:N-1
        NONL1(I,J)=1/2/DX*(-1/H(I+1,J)/2/DY*ZET(I+1,J)*(PSI(I+1,J+1)-PSI(I+1,J-1))+1/H(I-1,J)/2/DY*ZET(I-1,J)*(PSI(I-1,J+1)-PSI(I-1,J-1)));
        NONL(I,J)=NONL1(I,J)+1/2/DY*(1/H(I,J+1)/2/DX*ZET(I,J+1)*(PSI(I+1,J+1)-PSI(I-1,J+1))-1/H(I,J-1)/2/DX*ZET(I,J-1)*(PSI(I+1,J-1)-PSI(I-1,J-1)));
        BETA(I,J)=BETA(J)/H(I,J)*(PSI(I+1,J)-PSI(I-1,J))/2/DX;
        CURL(I,J)=1/2/DX/RHO*(TAUY(I+1,J)/H(I+1,J)-TAUY(I-1,J)/H(I-1,J))-1/2/DY/RHO*(TAUX(I,J+1)/H(I,J+1)-TAUX(I,J-1)/H(I,J-1));  
        DISS(I,J)= AH/DX^2*(ZET(I+1,J)-2*ZET(I,J)+ZET(I-1,J))+AH*(ZET(I,J+1)-2*ZET(I,J)+ZET(I,J-1))/DY^2;
        DIV1(I,J)=F(J)/2/DX*(-1/H(I+1,J)/2/DY*(PSI(I+1,J+1)-PSI(I+1,J-1))+1/H(I-1,J)/2/DY*(PSI(I-1,J+1)-PSI(I-1,J-1)));
        DIV(I,J)=DIV1(I,J)+F(J)/2/DY*(-1/H(I,J+1)/2/DX*(PSI(I+1,J+1)-PSI(I-1,J+1))-1/H(I,J-1)/2/DX*(PSI(I+1,J-1)-PSI(I-1,J-1))); 
            %Fixed Error - forgot "1/" for "1/H(I,J+1)"
            %Fixed Error - should subtract second half of equation
        ZN(I,J)=-NONL(I,J)-DIV(I,J)-BETA(I,J)+CURL(I,J)+DISS(I,J);    
        end
        end
        
        %%%%%%%% USE SOR TO TRANSFORM ZN INTO PHIN   
        
        c1 = 2/DX^2/(H(I+1,J)+H(I,J));
        c2 = 2/DY^2/(H(I,J+1)+H(I,J));
        c3 = 2/DX^2/(H(I,J)+H(I-1,J));
        c4 = 2/DY^2/(H(I,J)+H(I,J-1));
        c0 = c1 + c2 + c3 + c4;
        
        tic
        while anorm>anormf*eps_sor
        anorm=0;
            for i=2:M-1
            for j=2:N-1
                resid = c1*PHIN(i+1,j)+c2*PHIN(i,j+1)+c3*PHIN(i-1,j)+c4*PHIN(i,j-1)-c0*PHIN(i,j)-ZN(i,j);
                PHIN(i,j)=PHIN(i,j)+OMEGA/c0*resid;
                anorm=anorm+abs(resid);
            end;
            end;
            if anormf==0 anormf=anorm; end;
            itt=itt+1;
            fprintf('anorm=%e , anormf=%e itt=%i\n',anorm,anormf*eps_sor,itt);
%             if rem(itt,300) == 0
%                 figure(1)
%                 contourf(PHIN)
%                 colorbar
%                 title(['PHIN itt=' num2str(itt) ' itime=' num2str(itime)])
%                 xlabel('x')
%                 ylabel('y')
%             end    
        end
        toc

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for I=2:M-1
        for J=2:N-1
        PSINEW(I,J)=PSI(I,J)+DT*PHIN(I,J);
        ZETNEW(I,J)=ZET(I,J)+DT*ZN(I,J);
        end
        end  
        
        PSIOLD=PSI;
        ZETOLD=ZET;
        PSI=PSINEW;
        ZET=ZETNEW;
    
    else
        anormf=0.;
        anorm=1;
    
        for I=2:M-1
        for J=2:N-1
        NONL1(I,J)=1/2/DX*(-1/H(I+1,J)/2/DY*ZET(I+1,J)*(PSI(I+1,J+1)-PSI(I+1,J-1))+1/H(I-1,J)/2/DY*ZET(I-1,J)*(PSI(I-1,J+1)-PSI(I-1,J-1)));
        NONL(I,J)=NONL1(I,J)+1/2/DY*(1/H(I,J+1)/2/DX*ZET(I,J+1)*(PSI(I+1,J+1)-PSI(I-1,J+1))-1/H(I,J-1)/2/DX*ZET(I,J-1)*(PSI(I+1,J-1)-PSI(I-1,J-1)));
        BETA(I,J)=BETA(J)/H(I,J)*(PSI(I+1,J)-PSI(I-1,J))/2/DX;
        CURL(I,J)=1/2/DX/RHO*(TAUY(I+1,J)/H(I+1,J)-TAUY(I-1,J)/H(I-1,J))-1/2/DY/RHO*(TAUX(I,J+1)/H(I,J+1)-TAUX(I,J-1)/H(I,J-1));  
            % Fixed Error - TAUX & TAUY switch positions
        DISS(I,J)= AH/DX^2*(ZET(I+1,J)-2*ZET(I,J)+ZET(I-1,J))+AH*(ZET(I,J+1)-2*ZET(I,J)+ZET(I,J-1))/DY^2;
        DIV1(I,J)=F(J)/2/DX*(-1/H(I+1,J)/2/DY*(PSI(I+1,J+1)-PSI(I+1,J-1))+1/H(I-1,J)/2/DY*(PSI(I-1,J+1)-PSI(I-1,J-1)));
        DIV(I,J)=DIV1(I,J)+F(J)/2/DY*(-1/H(I,J+1)/2/DX*(PSI(I+1,J+1)-PSI(I-1,J+1))-1/H(I,J-1)/2/DX*(PSI(I+1,J-1)-PSI(I-1,J-1))); 
            %Fixed Error - forgot "1/" for "1/H(I,J+1)"
            %Fixed Error - should subtract second half of equation
        ZN(I,J)=-NONL(I,J)-DIV(I,J)-BETA(I,J)+CURL(I,J)+DISS(I,J); 
        end
        end
        
        %%%%%%%% USE SOR TO TRANSFORM ZN INTO PHIN   
        
        c1=2/DX^2/(H(I+1,J)+H(I,J));
        c2=2/DY^2/(H(I,J+1)+H(I,J));
        c3=2/DX^2/(H(I,J)+H(I-1,J));
        c4=2/DY^2/(H(I,J)+H(I,J-1));
        c0 = c1 + c2 + c3 + c4;
        
        tic
        while anorm>anormf*eps_sor
        anorm=0;
            for i=2:M-1
            for j=2:N-1
                resid = c1*PHIN(i+1,j)+c2*PHIN(i,j+1)+c3*PHIN(i-1,j)+c4*PHIN(i,j-1)-c0*PHIN(i,j)-ZN(i,j);
                PHIN(i,j)=PHIN(i,j)+OMEGA/c0*resid;
                anorm=anorm+abs(resid);
            end;
            end;
            if anormf==0 anormf=anorm; end;
            itt=itt+1;
            fprintf('anorm=%e , anormf=%e itt=%i\n',anorm,anormf*eps_sor,itt);
%             if rem(itt,300) == 0
%                 figure(1)
%                 contourf(PHIN)
%                 colorbar
%                 title(['PHIN itt=' num2str(itt) ' itime=' num2str(itime)])                 
%                 xlabel('x')
%                 ylabel('y')
%             end  
        end
        toc
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for I=2:M-1
        for J=2:N-1
            
        PSINEW(I,J)=PSIOLD(I,J)+2*DT*PHIN(I,J);  
        ZETNEW(I,J)=ZETOLD(I,J)+2*DT*ZN(I,J);
   
        end
        end
           
        PSIOLD=PSI;
        ZETOLD=ZET;
        PSI=PSINEW;
        ZET=ZETNEW;
   
    end

figure('Name','Barotropic Model','NumberTitle','off')
pcolor(LON,LAT,PSI)
title('PSI')
xlabel('X')
ylabel('Y')
shading interp
colorbar
    
filename = ['/media/473975BE618B0748/barotropic_model/Psi_timestep_' num2str(itt)];
eval(['print -djpeg ' filename])
clear filename

end




