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

data_bedprofile = load('prof1018.txt'); 


% constants
g = 9.81;                    % acceleration of gravity (m/s^2)
rho = 1025;                  % water density (kg/m^3)
fs = 2;                      % sampling frequency (Hz)
Npos = 5;                    % number of cross-shore positions considered
duration = length(data_lowtide)/fs; %seconds

length_data = length(data_lowtide(:,1)); %Length data set 
%Creating time vector 
time = linspace(0.25,duration,length_data);


% Initialization vectors for low tide 
% These vectors will be used to store the wave statistics at each position
Hrms_tot = zeros(Npos,1);  % root mean square height (m)
H13_tot  = zeros(Npos,1);  % significant wave height (m)
Hm_tot   = zeros(Npos,1);  % mean wave height (m)

% --------------------------------------
%     Computation of wave statistics
% --------------------------------------

%Loop for low tide only 

for ii=1:Npos  % loop on the positions
    
    low_tide_individual = zero_crossing(data_lowtide(:,ii),fs);
    
    Hm_tot(ii) = mean(low_tide_individual(:,1)); 
    H13_tot(ii) = significant_height(low_tide_individual(:,1)); 
    Hrms_tot(ii) = rms_height(low_tide_individual(:,1)); 
    
end

% --------------------------------------
%                  Output
% --------------------------------------

% Visualization of mean wave height, significant wave height and rms wave
% height for different positions cross shore of low tide 
figure() 
subplot(2,1,1)
plot(Hm_tot,'-*');
hold on 
plot(H13_tot,'-*'); 
plot(Hrms_tot,'-*'); 
legend('Mean wave height', 'Significant wave height', 'RMS wave height')
title('Mean height, significant wave height and rms wave height as a function of cross-shore position');
ylabel('h [m]'); 
grid on 
text(x,y,labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
%xlabel('Cross-shore Position'); 
%xticks([1 2 3 4 5])
%xticklabels({'P1','P3','P4','P5','P6'}); 
subplot(2,1,2)
plot(data_bedprofile(:,2)); 
xlabel('Cross-shore Position');
xlim([0 2520]); 
hold on 

%plot(data_bedprofile(:,1)); 



%% We compute for all cases (low, mid and high tides) 
clear ii 
clear jj 


% These vectors will be used to store the wave statistics at each position
Hrms_total = zeros(Npos,3);  % root mean square height (m)
H13_total  = zeros(Npos,3);  % significant wave height (m)
Hm_total   = zeros(Npos,3);  % mean wave height (m)

%We create an array with the names 
array_names = ["lowTide.txt", "midTide.txt", "highTide.txt"]; 
n = length(array_names); %Number of datasets 

%We run for all datasets 
for ii=1:Npos  % loop on the positions
    for jj= 1:n
        name = array_names(jj);
        wave_data = load(name); 
        waves_individual = zero_crossing(wave_data(:,ii),fs);
    
        Hm_total(ii,jj) = mean(waves_individual(:,1)); 
        H13_total(ii,jj) = significant_height(waves_individual(:,1)); 
        Hrms_total(ii,jj) = rms_height(waves_individual(:,1)); 
        
        
               
    end
    
end


%Plot of significant wave height as a function of root mean square wave
%height 

%We do polyfit of all dataset 

fit= polyfit(Hrms_total,H13_total,1);

%Theoretical relation for 1/3 and rms
theoretical_fit = sqrt(2)*Hrms_total; 

figure()
plot(fit);
hold on 
plot(theoretical_fit)

%Theoretical relation for mean and rms

Hm_theoretical = 0.89*Hrms_total; 

%Real 
H_div = Hm_total/Hrms_total; 










% ?