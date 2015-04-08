% heat equation problem
% solve T,t = 3T,xx + 4(x+1)(T-30)
% IC: T(x,0) = 25
% BC1: T(0,t) = 25 + 0.3t
% BC2: T,x(1,t) = 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define key variables
N = 10; % number of points on bar
dt = 0.00205; % time step - unstable if greater that roughly .00205 for this example
dx = 1/(N-1); % incremental length along bar of unit length
lambda = dt*3/dx^2; % k/C = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Told = zeros(10,1);
Tnew = zeros(10,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize T
for i = 1:N
	Told(i) = 25; % you could also do Told=25 for this problem
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0; %this is the initial time
tmax = 3; % this is the simulation time
% k =1;
fig1 = plot(Told);
axis([1,10,0,30])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while t < tmax
	set(fig1,'YData',Told);
	drawnow
	%pause(0.25)
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	% now start the forward Euler timestepping	
	t = t + dt; % update the time
	Tnew(1) = 25 + 0.3*t; % this is the LHS boundary condition
	for i = 2:N-1
		xi = (i-1)*dx; % this is the location of the ith node
		Tnew(i) = Told(i)+lambda*(Told(i+1)+Told(i-1)-2*Told(i))+4*dt*(xi+1)*(Told(i)-30);
	end
	TNp1 = 2*dx*4+Told(N-1); % this is the derivative condition at the end
	Tnew(N) = Told(N)+lambda*(TNp1+Told(N-1)-2*Told(N))+4*dt*((N-1)*dx+1)*(Told(i)-30);
	Told = Tnew;
end
