% One dimension model for vertical diffusion of heat and salt in a lake

clear all
close all

full_output = 1;%...................Toggle plots on/off  (1 = on, otherwise off) 

%% Parameters

% Physical

H = 250 ;%..........................depth of lake [m]
alphaT = 2e-4 ;%....................coefficient of thermal expansion [1/C]
alphaS = 7.6e-4 ;%..................coefficient of haline contraction [1/ppt]
t50 = 26.3 ;%.......................temperature @ 50m depth [C] <-
s50 = 32 ;%.........................salinity @ 50m depth [PSU]  <- From CTD profiles
rho50 = 1020.75 ;%..................density @ 50m depth [kg/m3] <-
lambda = 2.5e6 ;%...................latent heat of vaporization of water [J/kg]
beta = .08 ;%.......................Bowen ratio of sensible heat flux to latent heat flux
dQdT = 40 ;%........................change in heat flux with air-sea T difference [W/(m2*T)]
Cp = 4000 ;%........................specific heat capacity [J/(kg*degC]
kappa = 3.0e-6 ;%...................thermal diffusivity constant [m2/s] 
kappaS = kappa/100;%................salt diffusivity constant [m2/s]
g = 9.82;%..........................Graviational acceleration [m/s2]

% Numerical

nz = 250;%.........................# of grid points in z-direction
dz = H/nz;%........................vertical grid spacing
%dt = dz^2/(4*kappa);%..............Compute stable time step
dt = 3600;%........................Chosen time step ~ 1hr
years = 4;%........................Number of years to run model
nt = years*365*86400/dt + 5000;%..........# of time steps to compute
time=[1:nt]*dt;%...................Create time vector


%% Intitialize Matrices

T = zeros(nz,nt) + t50;%.....................Initialize temperature matrix
S = zeros(nz,nt) + s50;%.....................Initialize salinity matrix
Rho = zeros(nz,nt) + rho50;%.................Initialize density matrix
K = zeros(nz,nt) + kappa;%...................Initialize thermal diffusivity matrix
N2 = zeros(nz,nt);%..........................Initialize bouyancy frequency matrix

%% Boundary Conditions

%%% . '  , * . '',%%%
%%%   Salt Flux   %%%
%%% ' ' ,. ' * ', %%%

% Evaporation - using energy balance method -> Net Radiation = Sensible Heat Flux + Latent Heat Flux + Energy Storage Term
%               (Vallet-Couloumb, 2001)     -> Assume no change in energy storage (i.e. Energy Storage Term = 0)

day=[1:365*years];
Rmean = 250;%......................................Yearly average net insolation at Makapulapai Station [W/m2]
Ramp = 20;%........................................Amplitude of the yearly cycle at Makapulapai Station [W/m2]
Rnet = Ramp*sin(2*pi*day/365) + Rmean;%............Create yearly cycle [W/m2]
E = Rnet./(lambda*(beta+1)).*86400;%...............Calculate evaporation cycle [mm/day]
Expansion = round(nt/length(E))+1;%................Determine expansion factor needed for kron
E = kron(E,ones(1,Expansion)).*(1/Expansion);%.....Expand E so length(E)=nt, then scale down
E(nt+1:end) = [];%.................................Delete extra terms from expansion

%%% AXEL MOD
%E=5*(E-mean(E));

% Precipitation

%P = xlsread('/home/mike/modeling/lake_model_project/kauhakoprofilesandweatherstation/2011_P.xls', 1, 'basic')';%...................2011
%P = xlsread('/home/mike/modeling/lake_model_project/kauhakoprofilesandweatherstation/8-2002-8-2003_drought_P.xls', 1, 'basic');%...Drought
%P = xlsread('/home/mike/modeling/lake_model_project/kauhakoprofilesandweatherstation/8-2003-8-2004_rainy_P.xls', 1, 'basic')';%....Rainy
%Pinput = xlsread('/home/mike/modeling/lake_model_project/kauhakoprofilesandweatherstation/2008-2011-P.xls', 1, 'basic')';...........2008-2012
Pinput = xlsread('/home/akule/grissom7/Desktop/2008-2011-P.xls', 1, 'basic')';...........2008-2012
vi = find(isnan(Pinput));%...........................Find missing data
Pinput(vi) = 0;%.....................................Fill with zero precipitation
%P = repmat(P,1,years)%..............................Repeat 1 yr rainfall pattern for a number of years years
Expansion = round(nt/length(Pinput))+1;%.............Determine expansion factor needed for kron (data points per day)
P = kron(Pinput,ones(1,Expansion)).*(1/Expansion);%..Expand P so length(P)=nt, then scale down
P(nt+1:end) = [];%...................................Delete extra terms from expansion so length(P)=nt

%P = 0;

FWflux = E - P;%.....................................Calculate fresh water flux out [mm/hour]

 
Sflux = FWflux *1e-3/3600 *17/2;%....................Salinity flux = S * (E - P) * (mm -> m conversion) / mixed layer depth [PSU/hour]
Sflux = Sflux -mean(Sflux);


S(nz,2:end)=18 + Sflux(1,2:end);%....................Upper Boundary Condition - Prescribed SSS for all time
S(1,:)=s50;%.........................................Lower Boundary Condition
S(nz-4:nz,1)=[32 24 20 18 18];%......................Initial vertical profile S (from CTD)

% figure('Name','Salinity Boundary Conditions','NumberTitle','off')
% subplot(2,1,1)
% plot(time/86400/365,E-P)
% xlabel('Year')
% ylabel('mm/day')
% title('Evaporation - Precipitation')
% subplot(2,1,2)
% plot(time/86400/365,Sflux)
% title('Salinity flux')
% xlabel('Year')
% ylabel('PSU/hr')

%%% '''''''''''''''''%%%
%%% Temperature Flux %%%
%%% '''''''''''''''''%%%

% Import average daily temperature

%Tair = xlsread('/home/mike/modeling/lake_model_project/kauhakoprofilesandweatherstation/DailyT-2011.xls', 1, 'basic')';%.............2011
%Tair = xlsread('/home/mike/modeling/lake_model_project/kauhakoprofilesandweatherstation/DailyT-8-2002-8-2003.xls', 1, 'basic')';%...Drought
%Tair = xlsread('/home/mike/modeling/lake_model_project/kauhakoprofilesandweatherstation/DailyT-8-2003-8-2004.xls', 1, 'basic')';%...Rainy
%Tairinput = xlsread('/home/mike/modeling/lake_model_project/kauhakoprofilesandweatherstation/2008-2011-Tave.xls', 1, 'basic')';%...2008-2011
Tairinput = xlsread('/home/akule/grissom7/Desktop/2008-2011-Tave.xls', 1, 'basic')';%...2008-2011
vi = find(isnan(Tairinput));%.......................Find missing data
Tairinput(vi) = nanmean(Tairinput);%................Fill with mean temperature
%Tair = repmat(Tair,1,years);%......................Repeat 1 yr rainfall pattern for a number of years
Expansion = round(nt/length(Tairinput))+1;%.........Determine expansion factor needed for kron (data points per day)
Tair = kron(Tairinput,ones(1,Expansion));%..........Expand Tair so length(Tair)=nt
Tair(nt+1:end) = [];%...............................Delete extra terms from expansion


T(nz,:) = mean(Tair);%...................................Upper Boundary Condition - Prescribed SST for all time
T(1,2:end)=t50;%...................................Lower Boundary Condition T
T(nz-6:nz,1)=[26.7 27 28.5 33 29.5 28 28];%........Initial vertical profile T (from CTD)

% plot(time/86400/365,Tair)
% title('2008-2012 daily average temperature')
% xlabel('Year')
% ylabel('DegC')
   
%%%------------%%%
%%%   Density  %%%
%%%------------%%%

Rho(nz,2:end) = rho50*(ones(1,nt-1) + alphaS.*(S(nz,2:end)-s50) - alphaT*(T(nz,2:end) - t50));%....Upper Boundary Condition
Rho(nz-4:nz,1)=[1020 1014 1011 1009 1009];%........................................................Initial vertical profile Rho (from CTD)

for z = 2:nz-1%....................................................................................Initial thermal diffusivity profile
K(z,1) = kappa * (1 + 0.5*tanh(-(Rho(z+1,1) - Rho(z-1,1))/2/dz));
end

% plot(flipud(K),[-1:-1:-250])
% xlabel('Diffusivity [m2/s]')
% ylabel('Depth [m]')
% title('Initial Thermal Diffusivity Profile')


    for t=1:nt-1 %.............. Loop over all time except last time step
        
        T(250,t+1) = T(250,t) - dt* (T(250,t)-Tair(t))/(86400*3);
        S(250,t+1) = S(250,t) + dt* Sflux(t);
        Rho(250,t+1) = rho50 * (1 + alphaS*(S(250,t+1)-s50) - alphaT*(T(250,t+1) - t50));
        
        for z=2:nz-1 %...........Forward in Time, Centered in Space Scheme
            
            
            T(z,t+1) = T(z,t) + (dt/4/(dz^2)) * (K(z+1,t) - K(z-1,t)) * (T(z+1) - T(z-1)) + (K(z,t)*dt/(dz^2)) * (T(z+1,t) - 2*T(z,t) + T(z-1,t)) ;
            
            S(z,t+1) = S(z,t) + (kappaS*dt/(dz^2))*(S(z+1,t) - 2*S(z,t) + S(z-1,t));
            
            Rho(z,t+1) = rho50 * (1 + alphaS*(S(z,t+1)-s50) - alphaT*(T(z,t+1) - t50));
            
            K(z,t+1) = kappa * (1 + tanh(-(Rho(z+1,t) - Rho(z-1,t))/2/dz));
            
            N2(z,t+1) = -g / rho50 * (Rho(z+1,t) - Rho(z-1,t))/2/dz;
            
            dt = 0.25/max(K(:,t+1));
            
            time(t+1) = time(t) + dt;
            
        end
    end
    
    dt_stable = 0.5 / max(max(K))
    dt
    
%%% Various plots & movies

% for i = 1:nt
% plot(T(201:250,100*i),-50:-1)
% axis([20 34 -50 0])
% xlabel('Temperature [DegC]')
% ylabel('Depth [m]')
% title('Temperature Profile')
% pause(0.2)                 
% end
% 
% 
% for i = 1:nt
% plot(S(201:250,100*i),-50:-1)
% xlabel('Salinity [PSU]')
% ylabel('Depth [m]')
% title('Salinity Profile')
% pause(0.1)                 
% end
% 
% 
% for i = 1:nt
% plot(Rho(201:250,100*i),-50:-1)
% xlabel('Density [kg/m3]')
% ylabel('Depth [m]')
% title('Density Profile')
% pause(0.1)                 
% end

if full_output == 1

figure('Name','Evaporation','NumberTitle','off')
plot(time/86400/365,E)
title(['Evaporation over lake surface ~' num2str(round(sum(E)/years)) 'mm/yr'])
xlabel('Year')
ylabel('mm/hr')

figure('Name','Hoffmuller Diagrams - T,S,Rho','NumberTitle','off')
subplot(3,1,1)
pcolor(time/(86400*365),[-50:-1],T(201:250,:))
shading interp
colorbar
ylabel('Depth [m]')
title('T [degC]')

subplot(3,1,2)
pcolor(time/86400/365,[-50:-1],S(201:250,:))
shading interp
colorbar
ylabel('Depth [m]')
title('S [psu]')

subplot(3,1,3)
pcolor(time/(86400*365),[-50:-1],Rho(201:250,:))
hold on
shading interp
colorbar
hold on

xlabel('Year')
ylabel('Depth [m]')
title('Rho [kg/m3]')

figure('Name','Hoffmuller Diagrams - N2,K','NumberTitle','off')
subplot(2,1,1)
pcolor(time/(86400*365),[-50:-1],K(201:250,:))
shading interp
colorbar
ylabel('Depth [m]')
title('K [m2/s]')

subplot(2,1,2)
pcolor(time/(86400*365),[-50:-1],N2(201:250,:))
shading interp
colorbar
ylabel('Depth [m]')
xlabel('Time')
title('N2 [1/s]')

 
% 
%  figure('Name','Flux Terms','NumberTitle','off')
%  subplot(3,1,1)
%  plot(time/(86400*365),Tair);
%  ylabel('DegC')
%  title('Temperature at Makapulapai')
% % 
%  subplot(3,1,2)
%  plot(time/(86400*365),T(nz,:));
%  ylabel('DegC')
%  title('Prescribed Surface Temperature')
% % 
%  subplot(3,1,3)
%  plot(time/(86400*365),dt * Tflux)
%  xlabel('Year')
%  ylabel('DegC')
%  title('Temperature Change due to Heat Flux')
% 
% 
 figure('Name','Intial T,S,Rho Profiles','NumberTitle','off')
 subplot(1,3,1)
 plot(flipud(T(:,1)),[-1:-1:-nz])
 hold on
 plot(flipud(mean(T,2)),[-1:-1:-nz],'r')
 xlabel('Temperature [DegC]')
 ylabel('Depth [m]')
 title('Temperature Profile')
% 
 subplot(1,3,2)
 plot(flipud(S(:,1)),[-1:-1:-nz])
 hold on
 plot(flipud(mean(S,2)),[-1:-1:-nz],'r')
 xlabel('Salinity [PSU]')
 ylabel('Depth [m]')
 title('Salinity Profile')
% 
 subplot(1,3,3)
 plot(flipud(Rho(:,1)),[-1:-1:-nz])
 hold on
 plot(flipud(mean(Rho,2)),[-1:-1:-nz],'r')
 xlabel('Density [kg/m3]')
 ylabel('Depth [m]')
 title('Density Profile')
%  
% 
 figure('Name','Surface Values','NumberTitle','off')
 subplot(3,1,1)
 plot(time/(86400*365),T(nz,:))
 ylabel('Temperatture')
 title('Prescribed Surface Values')
% 
 subplot(3,1,2)
 plot(time/(86400*365),S(nz,:))
 ylabel('Salinity')
 
 subplot(3,1,3)
 plot(time/(86400*365),Rho(nz,:))
 xlabel('Year')
 ylabel('Density')

end