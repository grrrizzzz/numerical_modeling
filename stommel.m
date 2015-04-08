n=100;
dx=1/n;
a=zeros(n+1,n+1);
a(1,1)=1;
a(n+1,n+1)=1;
b=ones(n+1,1)*(-1);
b(1)=0;
b(n+1)=0;
for j=1:3
	eps=100*10^-j;
	for i=2:n
	a(i,i-1)=eps/dx^2-1/(2*dx);
	a(i,i+1)=eps/dx^2+1/(2*dx);
	a(i,i)=-2*eps/dx^2;
	end
	psi=inv(a)*b;
	plot(psi)
	pause(0.5)
end

