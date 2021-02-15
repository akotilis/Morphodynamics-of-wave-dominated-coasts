%% Section 1.2.2 -- Wave statistics 

%%  Wave statistics at low tide 

% Wave-by-wave analysis of the low tide data, Egmond Coast 3d dataset
% The objective of the script is to compute wave statistics at different
% locations along a cross-shore transect and examine their cross-shore
% evolution


% -------------------------------------
%            Initialisation
% -------------------------------------
clear all
close all

% load data
data_lowtide = load('lowTide.txt');
data_hightide = load('highTide.txt'); 
data_midtide = load('midTide.txt'); 


% constants
g = 9.81;                    % acceleration of gravity (m/s^2)
rho = 1025;                  % water density (kg/m^3)
fs = 2;                      % sampling frequency (Hz)
Npos = 5;                    % number of cross-shore positions considered
duration = length(data_lowtide)/fs; %seconds
length_data = length(data_lowtide(:,1)); 
time = linspace(0.25,duration,length_data);
% More?

% Initialisation vectors 
% These vectors will be used to store the wave statistics at each position
Hrms_tot = zeros(Npos,1);  % root mean square height (m)
H13_tot  = zeros(Npos,1);  % significant wave height (m)
Hm_tot   = zeros(Npos,1);  % mean wave height (m)

% --------------------------------------
%     Computation of wave statistics
% --------------------------------------


for ii=1:Npos  % loop on the positions
    
    low_tide_individual = zero_crossing(data_lowtide(:,ii),fs);
    
    Hm_tot(ii) = mean(low_tide_individual(:,1)); 
    H13_tot(ii) = significant_height(low_tide_individual(:,1)); 
    Hrms_tot(ii) = rms_height(low_tide_individual(:,1)); 
    
end

% --------------------------------------
%                  Output
% --------------------------------------

% visualisation of outputs
figure() 
plot(Hm_tot,'*');
hold on 
plot(H13_tot,'*'); 
plot(Hrms_tot,'*'); 
legend('Mean wave height', 'Significant wave height', 'RMS wave height')
title('Mean, significant wave height and rms wave height as a function of position'); 


%% We compute for all cases (low, mid and high tides) 

% These vectors will be used to store the wave statistics at each position
Hrms_total = zeros(Npos,3);  % root mean square height (m)
H13_total  = zeros(Npos,3);  % significant wave height (m)
Hm_total   = zeros(Npos,3);  % mean wave height (m)

array_names = ["lowTide.txt", "midTide.txt", "highTide.txt"]; 
n = length(array_names); %Number of datasets 

for ii=1:Npos  % loop on the positions
    for jj= 1:n
        name = array_names(jj);
        wave_data = load(name); 
        waves_individual = zero_crossing(wave_data(:,ii),fs);
    
        Hm_tot(ii,jj) = mean(low_tide_individual(:,1)); 
        H13_tot(ii,jj) = significant_height(low_tide_individual(:,1)); 
        Hrms_tot(ii,jj) = rms_height(low_tide_individual(:,1)); 
        
        fit_wv_rms = polyfit(H13_tot(:,jj),Hrms_tot(:,jj),2);
    end
    
end

%Plot of significant wave height as a function of root mean square wave
%height 

fit_wv_rms = polyfit(H13_tot,Hrms_tot,2);

%Theoretical relation for 1/3 and rms

H_theoretical = sqrt(2)*Hrms_tot; 


%Theoretical relation for mean and rms 

Hm_theoretical = 0.89*Hrms_tot; 








% ?