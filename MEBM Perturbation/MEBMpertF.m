clear all;
% Code for Moist/Dry Energy Balance Model

% This code solves the energy balance model used in Roe et al. (Nat. Geosci., 2015)
% The model operates in perturbation mode (i.e., the anomaly with respect to a prescribed climatology):
% You can specify:-
% the pattern of climate forcing, Rf
% the pattern of climate feedbacks, -B
% the ocean heat uptake, G
% the diffusivity, D
% whether you diffuse moist static energy, or just sensible heat

%time step in fraction of year
delt=1./500000; disp(['delt = ' num2str(delt)]) %1./50000;
NMAX=10000000; disp(['NMAX = ' num2str(NMAX)]) %60000

%set up x array (latitude).
jmx=101; %was 101
delx = 2.0/jmx;
x = [-1.0+delx/2:delx:1.0-delx/2]';x = x(:);
phi = asin(x)*180/pi;

% magnitude of diffusivity [units?]
Dmag = 0.2598; disp(['D = ' num2str(Dmag) ' W/(m2 K)'])% D = 0.2598 W/(m2 K) is the value used by TF10 
D=Dmag*ones(jmx+1,1); D=D(:); % diffusivity
 
% I think this C = rho * c * h_ml /(pi*1e7).
% this is consistent with a ~1m layer of dirt
% Note - heat capacity over LAND LAND LAND.
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
 
%% load in the climatological temperature from ERA-Interim for the control climate
load ERAtemperature.mat
T_ctrl = T;  % import observed T; units K
T_ctrl = 0.5*(T_ctrl+flipud(T_ctrl)); % average N & S hemispheres for symmetry
T_ctrl = interp1(sind(lat),T_ctrl,x,'linear');
q = eps*relhum/psfc*e0*exp(a*T_ctrl./(b+T_ctrl)); q=q(:);% here T is in oC. q is g kg-1
theta_e_ctrl = 1/cp*(cp*(T_ctrl+273.15) + L*q); % theta is mse divided by cp; note units of Kelvin are needed!!!
clear T lat;

%% --------------------------------------------------------------
% Choose a pattern of feedbacks (-B), forcing and ocean heat uptake.
% Default is CMIP5-mean pattern, but you can uniform feedback and forcing below

%CMIP5 ensemble-mean feedback and forcing values from 4xCO2 simulations (taken at year 100)
load CMIP5_Rf_G_lambda.mat %feedback, forcing, and heat uptake for 11 models
B = -interp1(CMIP5_lat,mean(CMIP5_lambda,2),phi,'linear');
R_frc = interp1(CMIP5_lat,mean((CMIP5_Rf + CMIP5_G),2),phi,'linear'); %CO2 forcing Rf plus ocean heat uptake G

figure(1);clf;
set(gcf,'defaultaxestickdir','out','defaultaxesticklength',[0.02 0],'defaultaxesbox','off')
set(gcf,'defaulttextfontname','times','defaultaxesfontname','times')
plot(CMIP5_lat,mean(CMIP5_lambda,2),'-k','linewidth',2)
xlim([-90 90])
ylim([-4 2])
set(gca,'fontsize',28)

figure(2);clf;
set(gcf,'defaultaxestickdir','out','defaultaxesticklength',[0.02 0],'defaultaxesbox','off')
set(gcf,'defaulttextfontname','times','defaultaxesfontname','times')
plot(CMIP5_lat,mean(CMIP5_Rf,2),'--b','linewidth',2)
hold on;
plot(CMIP5_lat,mean(CMIP5_G,2),'-g','linewidth',2)
xlim([-90 90])
ylim([-15 15])

%-------------------------------------------------
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Uncomment these to use flat forcing, uniform feedbacks, no ocean heat uptake
%R_frc = 7.8; % uniform forcing in [W/m2] for a quadrupling of CO2, taken as global average of CMIP5
%B = 1.4*ones(size(x)); B=B(:); disp('B and Rf are constant for expt.') % uniform feedback [W/m2/K], taken as average of CMIP5

%% Now solve climate change experiment
disp('**Doing climate change experiment**')
%--------------------------------------------------------------
%set up inital T profile 
T = 0.5*(1-1*x.^2); T = T(:);
Tinit=T;
Tglob = mean(Tinit);

% Timestepping loop
for n=1:NMAX
       
   
% spec. hum, and theta_e
   q = eps*relhum/psfc*e0*exp(a*(T+T_ctrl)./(b+(T+T_ctrl))); q=q(:);% here T is in oC. q is g kg-1
   theta_e = 1/cp*(cp*((T+T_ctrl)+273.15) + L*q); % note units of Kelvin are needed!!!
   theta_e_pert = theta_e-theta_e_ctrl; %anomalous theta_e

% Calculate new T.
% Diffuse anomalous moist static energy (theta_e_pert)
  dT = delt/Cl*(R_frc - (B.*T) + Mdiv*theta_e_pert); if n == 1, disp(['Diffusing theta_e anomaly']), end;

   T = T + dT; 
   
% Check to see if global mean energy budget has converged:   
   Fglob=mean(R_frc - (B.*T));
   if (abs(Fglob) < 0.001), break; end
end

divF_pert = -Mdiv*theta_e_pert;
h_pert = theta_e_pert*cp;

%% Plot some things

 xt = [-75 -60 -45 -30 -15 0 15 30 45 60 75];
 xt = sind(xt);
 xtl = ['-75'; '-60'; '-45'; '-30'; '-15'; '  0'; ' 15'; ' 30'; ' 45'; ' 60';' 75'];

 figure(1);clf;
plot(x,T,'linewidth',3); grid on
hold on;
plot(sind(CMIP5_lat),mean(CMIP5_T,2),'k','linewidth',3); %compare with CMIP5 T anomaly over same latitude ranges
xlabel('Latitude')
ylabel('\Delta T (^oC)')
set(gca,'XTick',xt);
set(gca,'XTicklabel',xtl);
set(gca,'fontsize',14,'XLim',[-1 1],'YLim',[0 15],'YTick',0:5:20)
 
figure(2);clf;
plot(x,divF_pert,'linewidth',3); grid on
xlabel('Latitude')
ylabel('\nabla . F (W m^{-2})')
set(gca,'XTick',xt);
set(gca,'XTicklabel',xtl);
set(gca,'fontsize',14)

figure(3);
plot(x,h_pert,'linewidth',3); grid on
xlabel('Latitude')
ylabel('\Delta h (J kg^{-1})')
set(gca,'XTick',xt);
set(gca,'XTicklabel',xtl);
set(gca,'fontsize',14)

 