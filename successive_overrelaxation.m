
clear

jmt=6;km=6;
Lx=0.5;
Ly=0.5;
dx=Lx/(jmt-1);dy=Ly/(km-1);
x=0:dx:Lx;
y=0:dy:Ly;
jmt=length(x);
km=length(y);

psi=zeros(jmt,km);
psi(:,end)=200*x';
psi(end,:)=200*y;




rja= (cos(pi/jmt)+(dx/dy)^2*cos(pi/km))/(1+(dx/dy)^2);
omega=2/(1+(1-rja^2)^0.5);

%%%%%%%%%%%% HERE IS where the interesting SOR part starts

cjp=1/dx^2 ;
cjm=1/dx^2 ;

ckp= 1/dy^2;
ckm= 1/dy^2;

e=-2/dx^2-2/dy^2 ;
f=0;

eps_sor=0.001;
anormf=0.;
anorm=1;

itt=0;

while anorm>anormf*eps_sor
 anorm=0;
 for j=2:jmt-1
  for k=2:km-1
    resid=psi(j+1,k)*cjp+psi(j-1,k)*cjm+psi(j,k+1)*ckp+psi(j,k-1)*ckm ...
           +e*psi(j,k) -f;
    psi(j,k)=psi(j,k)-omega*resid/e;
    anorm=anorm+abs(resid);
  end;
 end;
 if anormf==0 anormf=anorm; end;
 itt=itt+1;
 fprintf('anorm=%e , anormf=%e itt=%i\n',anorm,anormf*eps_sor,itt);

 clf
contourf(x,y,psi')
xlabel('x')
ylabel('y')
title('\nabla^2 \phi = 0')
colorbar 
pause(0.5)
 
end;

% clf
% contourf(x,y,psi')
% xlabel('x')
% ylabel('y')
% title('\nabla^2 \phi = 0')
% colorbar 