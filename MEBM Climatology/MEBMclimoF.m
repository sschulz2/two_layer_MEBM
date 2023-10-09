clear all;
% Code for Moist/Dry Energy Balance Model

% This code solves the energy balance model used in Roe et al. (Nat. Geosci., 2015)
% The model operates in climatology mode
% You can specify:-
% the insolation Q0
% the OLR parameters, A0,B0
% the diffusivity, D
% albedo of ocean and ice
% whether you diffuse moist static energy, or just sensible heat

%to compare to observed (ERA-I) T from HW1:
%load ERAtemperature.mat
%T_ERA = T;
%lat_ERA = lat;
%clear T lat;

%time step in fraction of year
delt=1./500000; disp(['delt = ' num2str(delt)]) %1./50000;
NMAX=10000000; disp(['NMAX = ' num2str(NMAX)]) %60000

%set up x array (latitude).
jmx=101; %was 101
delx = 2.0/jmx;
x = [-1.0+delx/2:delx:1.0-delx/2]';x = x(:);
phi = asin(x)*180/pi;

% climate parameters
Q0 = 342;                          % [W m-2]  solar constant n.b. 4Q0 = 1368
A0= 207.3;                         % Size of longwave cooling constant [W/m2]
B0 = 2.09;                         % [W m-2 degC-1] OLR response NOTE UNITS
alf_noice = 0.3;                   % [] ice free albedo.
alf_ice = 0.55; %0.55;                    % [] ice covered albedo.

A = A0;
B=B0*ones(size(x)); B=B(:);

% magnitude of diffusivity [units?]
Dmag = 0.2598; disp(['D = ' num2str(Dmag) ' W/(m2 K)'])% D = 0.2598 W/(m2 K) is the value used by TF10 
D=Dmag*ones(jmx+1,1); D=D(:); % diffusivity
 
% I think this C = rho * c * h_ml /(pi*1e7).
% this is consistent with a ~1m layer of dirt
% Note - heat capacity over LAND. -- low heat capacity to reach equilibrium quickly
Cl = 0.2; disp(['Cl = ' num2str(Cl) ' J /(m2 K)'])

% Moisture parameters
relhum = 0.8;   % relative humidity
eps = 0.622;    % moisture coonstant
psfc = 9.8e4;   % (Pa)
e0 = 611.2;     % vap. press (Pa)
a = 17.67; b = 243.5;   % sat vap constants !!T must be in temperature
L = 2.45e6;         % latent heat of vaporization (J kg-1)
cp = 1004;          % (J kg-1 K-1)

%-------------------------------------------------
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%For dry energy balance model, uncomment these values
% Dmag = 0.44; % magnitude of diffusivity [W/(m2 K)]
% D=Dmag*ones(jmx+1,1); D=D(:); % diffusivity
% relhum = 0;  % switch off humidity
% disp('Diffusing sensible heat only')
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%-------------------------------------------------
 
% Use setupfast to create a matrix that calculates D*d/dx[ (1-x^2)d/dx] of
% whatever it operates on.
 [~,Mdiv]=setupfastM(delx,jmx,D,0.,1.0,delt);
 
%% Now solve for climate
disp('**Solving for climate**')
%--------------------------------------------------------------
%set up inital T profile 
T = 0.5*(1-1*x.^2); T = T(:);
Tinit=T;
Tglob = mean(Tinit);

% Timestepping loop
for n=1:NMAX
       
% Calculate src for this loop.
% Source  (ASR)
alf = alf_noice*ones(size(x));idx = find(T<=-10); alf(idx) = alf_ice;
Src = Q0*(1-0.241*(3*x.^2-1)).*(1-alf);
   
% spec. hum, and theta_e
   q = eps*relhum/psfc*e0*exp(a*(T)./(b+(T))); q=q(:);% here T is in oC. q is g kg-1
   theta_e = 1/cp*(cp*((T)+273.15) + L*q); % note units of Kelvin are needed!!!

% Calculate new T.
% Diffuse moist static energy (theta_e)
  dT = delt/Cl*(Src -A - (B.*T) + Mdiv*theta_e); if n == 1, disp(['Diffusing theta_e anomaly']), end;

   T = T + dT; 
   
% Check to see if global mean energy budget has converged:   
   Fglob=mean(Src - A - (B.*T));
   if (abs(Fglob) < 0.001), break; end
end

Fglob

divF = -Mdiv*theta_e;
h = theta_e*cp;

%% Plot some things

 xt = [-75 -60 -45 -30 -15 0 15 30 45 60 75];
 xt = sind(xt);
 xtl = ['-75'; '-60'; '-45'; '-30'; '-15'; '  0'; ' 15'; ' 30'; ' 45'; ' 60';' 75'];

 figure(1);clf;
plot(x,T,'linewidth',3); grid on
hold on;
%plot(sind(lat_ERA),T_ERA,'linewidth',3);
xlabel('Latitude')
ylabel('\Delta T (^oC)')
set(gca,'XTick',xt);
set(gca,'XTicklabel',xtl);
set(gca,'fontsize',14)

figure(2);clf;
plot(x,divF,'linewidth',3); grid on
xlabel('Latitude')
ylabel('\nabla . F (W m^{-2})')
set(gca,'XTick',xt);
set(gca,'XTicklabel',xtl);
set(gca,'fontsize',14)

figure(3);
plot(x,h,'linewidth',3); grid on
xlabel('Latitude')
ylabel('\Delta h (J kg^{-1})')
set(gca,'XTick',xt);
set(gca,'XTicklabel',xtl);
set(gca,'fontsize',14)


Snk = A + B.*T;
figure(4)
plot(x,divF,'b-',x,Src,'g-',x,Snk,'r-',...
    x,Src-Snk-divF,'k','linewidth',1.5)
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl);
title('Terms in the energy balance','fontsize',12)
xlabel('Latitude','fontsize',14); ylabel('W m^{-2}','fontsize',14)
h1 = legend('Transport (\nabla F)','net SW','net OLR','Sum',...
    'location','west', 'orientation','vertical');
grid on

