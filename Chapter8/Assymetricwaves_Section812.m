%% 8.1.2 Assymetric waves 

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
PHI = 0;     % skewed wave 

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
title('Net sediment in the x-direction for assymetric waves (D50 = 0.3 mm)')
grid on 
xlabel('r')
ylabel('Q_s_,_x (m??2/s)') 
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


%% Now for fine sediment (D50 = 0.1 mm) 

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
title('Net sediment in the x-direction for assymetric waves (D50 = 0.1 mm)')
grid on 
xlabel('r')
ylabel('Q_s_,_x (m??2/s)')
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



