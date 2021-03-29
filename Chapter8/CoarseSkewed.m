% Computation of the sediment transport using the SANTOSS model for a
% coarse sand and skewed wave

% -------------------------------------------------
%                   INITIALISATION
% -------------------------------------------------
clear all
close all
clc


%%% Input parameters for the SANTOSS model

% Wave characteristics
Uw = 1.2;   % orbital velocity amplitude   [m/s]
T = 7;      % period [s]

% Wave shape parameters 
r = 0:0.15:0.6;  % define "r" vector from 0 to 0.6 with step size of 0.15
Nr = length(r);  % number of elements of r
PHI = -pi/2;     % skewed wave 

% Sediment characteristics
D50 = 0.3;  % D50 in mm
D90 = D50;  % D90 in mm
Rhos = 2650; % sediment density in kg/m^3

%%% Initialisation of the output vectors
Qsx = zeros(1,Nr);  % net sediment in the x-direction (m2/s)
Qsy = zeros(1,Nr);  % net sediment in the y-direction (m2/s)
Occ = zeros(1,Nr);  % dimensionless sed. load entrained during the crest period and transported during the crest
Oct = zeros(1,Nr);  % dimensionless sed. load entrained during the crest period and transported during the  trough
Otc = zeros(1,Nr);  % dimensionless sed. load entrained during the trough period and transported during the crest
Ott = zeros(1,Nr);  % dimensionless sed. load entrained during the trough period and transported during the trough

% -------------------------------------------------
%          COMPUTATION OF THE SED TRANSPORT
% -------------------------------------------------

for rI = 1:Nr     % loop on the different values of r considered
        
        % 1- computation of the time-series of orbital velocity 
         [u{rI} t{rI}] = waveshape(r(rI),PHI,Uw,T); 
        
        % 2- computation of the velocity skewness R and the acceleration skewness beta
         [R(rI) beta(rI) acceleration{rI}] = velocity_acceleration_skewness(u{:,rI},t{:,rI}); 

        
        % 3- computation of the root-mean squared orbital velocity Urms: should be in CM/S!
        %The root mean squared orbital velocity Urms can be computed as the standard 
        %deviation of the orbital velocity timeseries u(t).
        % We first change it from m/s to cm/s 
        u_cm{rI} = u{rI}.*100; 
        Urms(rI) = std(u_cm{rI}); 
        
        % 4- sediment transport calculation
        [Qsx(rI) Qsy(rI) Occ(rI) Oct(rI) Ott(rI) Otc(rI)] = SANTOSSmodel(D50,D90,Rhos,T,Urms(rI),R(rI),beta(rI),0,0);
end;


% -------------------------------------------------
%                 VISUALISATION
% -------------------------------------------------
figure() 
subplot(3,1,1) 
plot(r, Qsx, '*') 
title('Net sediment in the x-direction (D50 = 0.3 mm)')
grid on 
xlabel('r')
ylabel('Q_s_,_x (mˆ2/s)') 
subplot(3,1,2)
plot(r, Occ, '*k') 
title('Dimensionless sed. load entrained during the crest period and transported during the crest (Occ)')
grid on 
xlabel('r')
ylabel('Occ')
subplot(3,1,3)
plot(r,Ott, '*r') 
title('Dimensionless sed. load entrained during the trough period and transported during the trough (Ott)')
grid on 
xlabel('r')
ylabel('Ott')





%% Case with finer sediment 0.1 mm 


D50 = 0.1;  % D50 in mm


%%% Initialisation of the output vectors
Qsx = zeros(1,Nr);  % net sediment in the x-direction (m2/s)
Qsy = zeros(1,Nr);  % net sediment in the y-direction (m2/s)
Occ = zeros(1,Nr);  % dimensionless sed. load entrained during the crest period and transported during the crest
Oct = zeros(1,Nr);  % dimensionless sed. load entrained during the crest period and transported during the  trough
Otc = zeros(1,Nr);  % dimensionless sed. load entrained during the trough period and transported during the crest
Ott = zeros(1,Nr);  % dimensionless sed. load entrained during the trough period and transported during the trough

% -------------------------------------------------
%          COMPUTATION OF THE SED TRANSPORT
% -------------------------------------------------

for rI = 1:Nr     % loop on the different values of r considered
        
        % 1- computation of the time-series of orbital velocity 
         [u{rI} t{rI}] = waveshape(r(rI),PHI,Uw,T); 
        
        % 2- computation of the velocity skewness R and the acceleration skewness beta
         [R(rI) beta(rI) acceleration{rI}] = velocity_acceleration_skewness(u{:,rI},t{:,rI}); 

        
        % 3- computation of the root-mean squared orbital velocity Urms: should be in CM/S!
        %The root mean squared orbital velocity Urms can be computed as the standard 
        %deviation of the orbital velocity timeseries u(t).
        % We first change it from m/s to cm/s 
        u_cm{rI} = u{rI}.*100; 
        Urms(rI) = std(u_cm{rI}); 
        
        % 4- sediment transport calculation
        [Qsx(rI) Qsy(rI) Occ(rI) Oct(rI) Ott(rI) Otc(rI)] = SANTOSSmodel(D50,D90,Rhos,T,Urms(rI),R(rI),beta(rI),0,0);
end;

figure() 
subplot(5,1,1)
plot(r, Qsx, '*') 
grid on 
title('Net sediment in the x-direction (D50 = 0.1 mm)')
grid on 
xlabel('r')
ylabel('Q_s_,_x (mˆ2/s)')
subplot(5,1,2) 
plot(r, Occ, '*k') 
title('Dimensionless sed. load entrained during the crest period and transported during the crest (Occ)')
grid on 
xlabel('r')
ylabel('Occ')
subplot(5,1,3) 
plot(r,Oct, '*m')
title('Dimensionless sed. load entrained during the crest period and transported during the  trough (Oct)')
grid on 
xlabel('r')
ylabel('Oct')
subplot(5,1,4)
plot(r, Otc, '*b') 
title('Dimensionless sed. load entrained during the trough period and transported during the crest (Otc)')
grid on 
xlabel('r')
ylabel('Otc')
subplot(5,1,5)
plot(r, Ott, '*r')
title('Dimensionless sed. load entrained during the trough period and transported during the trough (Ott)')
grid on 
xlabel('r')
ylabel('Ott')


%% 8.2.1 Sediment transport analysis 

%Sediment transport for Egmond for midtide 
clear all; 


%%% Input parameters for the SANTOSS model

load('BattjesJansenmidtideChapter6.mat')

% Wave characteristics for midtide 
Uw = Uw(1:N_last);   % orbital velocity amplitude   [m/s]
T = 6.69;      % period [s]

% Wave shape parameters 
%r = 0:0.15:0.6;  % define "r" vector from 0 to 0.6 with step size of 0.15
Nr = length(r);  % number of elements of r
PHI = phi;     % skewed wave 

% Sediment characteristics
D50 = 0.225;  % D50 in mm for Egmond 
D90 = D50;  % D90 in mm
Rhos = 2650; % sediment density in kg/m^3

%%% Initialisation of the output vectors
Qsx = zeros(1,Nr);  % net sediment in the x-direction (m2/s)
Qsy = zeros(1,Nr);  % net sediment in the y-direction (m2/s)
Occ = zeros(1,Nr);  % dimensionless sed. load entrained during the crest period and transported during the crest
Oct = zeros(1,Nr);  % dimensionless sed. load entrained during the crest period and transported during the  trough
Otc = zeros(1,Nr);  % dimensionless sed. load entrained during the trough period and transported during the crest
Ott = zeros(1,Nr);  % dimensionless sed. load entrained during the trough period and transported during the trough

% -------------------------------------------------
%          COMPUTATION OF THE SED TRANSPORT
% -------------------------------------------------

for rI = 1:Nr     % loop on the different values of r considered
        
        % 1- computation of the time-series of orbital velocity 
         [u{rI} t{rI}] = waveshape(r(rI),PHI(rI),Uw(rI),T); 
         
        
        % 2- computation of the velocity skewness R and the acceleration skewness beta
         [R(rI) beta(rI) acceleration{rI}] = velocity_acceleration_skewness(u{:,rI},t{:,rI}); 

        
        % 3- computation of the root-mean squared orbital velocity Urms: should be in CM/S!
        %The root mean squared orbital velocity Urms can be computed as the standard 
        %deviation of the orbital velocity timeseries u(t).
        % We first change it from m/s to cm/s 
        u_cm{rI} = u{rI}.*100; 
        Urms(rI) = std(u_cm{rI}); 
        
        % 4- sediment transport calculation
        [Qsx(rI) Qsy(rI) Occ(rI) Oct(rI) Ott(rI) Otc(rI)] = SANTOSSmodel(D50,D90,Rhos,T,Urms(rI),R(rI),beta(rI),0,0);
end;

%Analyse the evolution of Qs,x as a function of x.

%Calculating dQs,x/dx 

dQsx = gradient(Qsx); 
dx = waves.x(2) - waves.x(1); 

dQsxdx = dQsx/dx; 

bed_profile = load('prof1018.txt'); 

%Figure of sediment transport, gradient of sediment transport and bed
%profile 
figure() 
subplot(3,1,1) 
plot(waves.x(1:N_last), Qsx) 
title('Net sediment transport in the x direction for Egmond (mid tide, D50 = 0.225 mm)') 
xlabel('x (m)')
ylabel('Q_s_,_x(mˆ2/s)')
grid on 
xlim([4000 max(waves.x(1:N_last))]) 
subplot(3,1,2) 
plot(waves.x(1:N_last),dQsxdx)
title('Gradient of sediment transport in the x direction for Egmond (mid tide D50 = 0.225 mm)')
xlabel('x (m)')
ylabel('dQ_x_,_s/dx (m/s)')
grid on 
xlim([4000 max(waves.x(1:N_last))]) 
subplot(3,1,3)
plot(waves.x,bed_profile(:,2),'k')
title('Bed profile') 
xlabel('x (m)')
ylabel('z (m)') 
grid on 
xlim([4000 max(waves.x(1:N_last))]) 

% Bed level change after 10 days 

epsilon = 0.6 %bed porosity 
deltat = 10*86400; %10 days in seconds (duration bed level change) 

num = -epsilon*deltat; 

deltazeta = dQsxdx*num; %bed level change (meters) 

figure() 
plot(waves.x,bed_profile(:,2)) 
hold on; 
plot(waves.x(1:N_last), deltazeta) 
grid on 
xlabel('x(m)') 
ylabel('z(m) and \delta z (m) ') 
legend('Bed profile', 'Bed level after 10 days')
xlim([4000 max(waves.x(1:N_last))]) 
ylim([-5 5])

%% 8.2.2 Sediment transport by waves and current 

%Initialization for calculation of the undertow 
E =  waves.E(1:N_last); %wave energy 
Er = waves.Er(1:N_last); %roller energy 
c = waves.c(1:N_last); %phase velocity 
h = waves.ht(1:N_last); %total water depth (m) 
rho = 1000; %density of water (kg/mˆ3) 
Hrms = waves.Hrms(1:N_last); %significant wave height (m) 


%%% Initialisation of the output vectors
Qsx_undertow = zeros(1,Nr);  % net sediment in the x-direction (m2/s)
Qsy_undertow = zeros(1,Nr);  % net sediment in the y-direction (m2/s)
Occ_undertow = zeros(1,Nr);  % dimensionless sed. load entrained during the crest period and transported during the crest
Oct_undertow = zeros(1,Nr);  % dimensionless sed. load entrained during the crest period and transported during the  trough
Otc_undertow = zeros(1,Nr);  % dimensionless sed. load entrained during the trough period and transported during the crest
Ott_undertow = zeros(1,Nr);  % dimensionless sed. load entrained during the trough period and transported during the trough

% -------------------------------------------------
%          COMPUTATION OF THE SED TRANSPORT
% -------------------------------------------------

for rI = 1:Nr     % loop on the different values of r considered
        
        %Computation of the undertow 
        Ux(rI) = magnitude_undertow(E(rI), Er(rI), c(rI), h(rI), rho, Hrms(rI));

        %We change from m/s to cm/s 
        Ux_cm(rI) = Ux(rI)*100; 
        
        % 4- sediment transport calculation
        [Qsx_undertow(rI) Qsy_undertow(rI) Occ_undertow(rI) Oct_undertow(rI) Ott_undertow(rI) Otc_undertow(rI)] = SANTOSSmodel(D50,D90,Rhos,T,Urms(rI),R(rI),beta(rI),0,Ux_cm(rI));
end;

figure()
subplot(5,1,1) 
plot(waves.x(1:N_last), Qsx)
hold on 
plot(waves.x(1:N_last), Qsx_undertow)
title('Net sediment transport in the x direction for Egmond (mid tide, D50 = 0.225 mm)') 
xlabel('x (m)')
ylabel('Q_s_,_x(mˆ2/s)')
grid on 
xlim([4000 max(waves.x(1:N_last))])
ylim([-2e-5 2e-5]) 
legend('Q_s_,_x waves', 'Q_s_,_x undertow')
%Occ
subplot(5,1,2) 
plot(waves.x(1:N_last),Occ) 
title('Dimensionless sed. load entrained during the crest period and transported during the crest (\Omega cc)')
grid on 
xlabel('x (m)')
ylabel('\Omega cc')
hold on
plot(waves.x(1:N_last),Occ_undertow) 
legend('\Omega cc waves', '\Omega cc undertow')
xlim([4000 max(waves.x(1:N_last))])
%Oct 
subplot(5,1,3) 
plot(waves.x(1:N_last),Oct) 
title('Dimensionless sed. load entrained during the crest period and transported during the  trough (\Omega ct)')
grid on 
xlabel('x (m)')
ylabel('\Omega ct')
hold on
plot(waves.x(1:N_last),Oct_undertow) 
legend('\Omega ct waves', '\Omega ct undertow')
xlim([4000 max(waves.x(1:N_last))])
ylim([-0.3 0.3]) 
%Otc 
subplot(5,1,4) 
plot(waves.x(1:N_last),Otc) 
title('Dimensionless sed. load entrained during the trough period and transported during the crest (\Omega tc)')
grid on 
xlabel('x (m)')
ylabel('\Omega tc')
hold on
plot(waves.x(1:N_last),Otc_undertow) 
legend('\Omega tc waves', '\Omega tc undertow')
xlim([4000 max(waves.x(1:N_last))])
%Ott 
subplot(5,1,5) 
plot(waves.x(1:N_last),Ott) 
title('Dimensionless sed. load entrained during the trough period and transported during the trough (\Omega tt)')
grid on
xlabel('x (m)')
ylabel('\Omega tt')
hold on
plot(waves.x(1:N_last),Ott_undertow) 
legend('\Omega tt waves', '\Omega tt undertow')
xlim([4000 max(waves.x(1:N_last))])
ylim([0 4]) 

%Gradient sediment transport undertow 

dQsx_undertow = gradient(Qsx_undertow); 
dx = waves.x(2) - waves.x(1); 

dQsxdx_undertow = dQsx_undertow/dx;

%Bed level after 10 days undertow 
deltazeta_undertow = dQsxdx_undertow*num; %bed level change (meters) 




%Add ticks to ylabel 
figure() 
subplot(2,1,1) 
plot(waves.x(1:N_last),dQsxdx)
hold on 
plot(waves.x(1:N_last),dQsxdx_undertow)
title('Gradient of sediment transport in the x direction for Egmond (mid tide D50 = 0.225 mm)')
xlabel('x (m)')
ylabel('dQ_x_,_s/dx (m/s)')
grid on 
xlim([4000 max(waves.x(1:N_last))]) 
legend('dQ_x_,_s/dx waves', 'dQ_x_,_s/dx undertow')
ylim([-1.2e-6 1.2e-6]) 
subplot(2,1,2)
plot(waves.x,bed_profile(:,2)) 
hold on; 
plot(waves.x(1:N_last), deltazeta) 
plot(waves.x(1:N_last), deltazeta_undertow) 
grid on 
xlabel('x(m)') 
ylabel('z(m) and \delta z (m) ') 
legend('Bed profile', 'Bed level after 10 days waves', 'Bed level after 10 days undertow')
xlim([4000 max(waves.x(1:N_last))]) 
ylim([-5 5])



%% 8.2.3 Influence of the offshore wave conditions

%% 8.2.2 Sediment transport by waves and current 

%Initialization for calculation of the undertow 
E =  waves.E(1:N_last); %wave energy 
Er = waves.Er(1:N_last); %roller energy 
c = waves.c(1:N_last); %phase velocity 
h = waves.ht(1:N_last); %total water depth (m) 
rho = 1000; %density of water (kg/mˆ3) 
Hrms = waves.Hrms(1:N_last); %significant wave height (m) 


%%% Initialisation of the output vectors
Qsx_undertow = zeros(1,Nr);  % net sediment in the x-direction (m2/s)
Qsy_undertow = zeros(1,Nr);  % net sediment in the y-direction (m2/s)
Occ_undertow = zeros(1,Nr);  % dimensionless sed. load entrained during the crest period and transported during the crest
Oct_undertow = zeros(1,Nr);  % dimensionless sed. load entrained during the crest period and transported during the  trough
Otc_undertow = zeros(1,Nr);  % dimensionless sed. load entrained during the trough period and transported during the crest
Ott_undertow = zeros(1,Nr);  % dimensionless sed. load entrained during the trough period and transported during the trough

% -------------------------------------------------
%          COMPUTATION OF THE SED TRANSPORT
% -------------------------------------------------

for rI = 1:Nr     % loop on the different values of r considered
        
        %Computation of the undertow 
        Ux(rI) = magnitude_undertow(E(rI), Er(rI), c(rI), h(rI), rho, Hrms(rI));

        %We change from m/s to cm/s 
        Ux_cm(rI) = Ux(rI)*100; 
        
        % 4- sediment transport calculation
        [Qsx_undertow(rI) Qsy_undertow(rI) Occ_undertow(rI) Oct_undertow(rI) Ott_undertow(rI) Otc_undertow(rI)] = SANTOSSmodel(D50,D90,Rhos,T,Urms(rI),R(rI),beta(rI),0,Ux_cm(rI));
end;


















