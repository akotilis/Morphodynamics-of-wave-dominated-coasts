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
bed_height = data_bedprofile(:,2);
bed_x = data_bedprofile(:,1);

% constants
g = 9.81;                    % acceleration of gravity (m/s^2)
rho = 1025;                  % water density (kg/m^3)
fs = 2;                      % sampling frequency (Hz)
Npos = 5;                    % number of cross-shore positions considered
duration = length(data_lowtide)/fs; %seconds
delta_t = 1/fs;
P_points = [4478 4765 4790 4814 4835];      %x-locations of the sensors
P_depths = bed_height(floor(P_points/2));   %depth of the sensors

%Length data set 
length_data = length(data_lowtide(:,1)); 

%Creating time vector 
time = linspace(delta_t,duration,length_data);


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
    H = low_tide_individual(:,1);
    Hm_tot(ii) = mean(H); 
    H13_tot(ii) = significant_height(H); 
    Hrms_tot(ii) = rms_height(H); 
    
end

% --------------------------------------
%                  Output
% --------------------------------------

% Visualization of mean wave height, significant wave height and rms wave
% height for different positions cross shore of low tide 
figure() 
subplot(2,1,1)
plot(P_points/1000,Hm_tot,'*');
hold on 
plot(P_points/1000,H13_tot,'*'); 
plot(P_points/1000,Hrms_tot,'*'); 
legend('H_{mean}', 'H_{1/3}', 'H_{rms}')
title('H_{mean}, H_{1/3} and H_{rms} as a function of cross-shore position');
xlabel('Cross-shore Position [km]'); 
ylabel('H [m]'); 
xlim([4.4 5]); 
grid on 
subplot(2,1,2)
plot(bed_x/1000,bed_height, 'b-');
hold on
plot(P_points/1000, P_depths, 'r.', 'linewidth',20);
txt = ['P1';'P3';'P4';'P5';'P6'];
text(P_points/1000,P_depths+1,txt)
xlabel('Cross-shore Position [km]');
ylabel('Bed height [m]')
title('Bed height and sensor positions')
xlim([4.4 5]); 
ylim([-10 2]);
grid on 



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


%% Plot of significant wave height as a function of root mean square wave height

%Using a reference height array for the fit
H_ref = [0 0.25 0.5 0.75 1 1.25 1.5 1.75 2]; 

%Using polyfit to find a linear relation between Hrms and H13
p = polyfit(Hrms_total, H13_total,1);

%Linear fit using the slope and intercept found by polyfit
slope = p(1);
intercept = p(2);
fit = slope*H_ref + intercept;

%Theoretical relation for H1/3 and Hrms
theoretical_fit = sqrt(2)*H_ref; 

%0.75,1.02
%1.25, 1.732
slope_theoretical = (1.732-1.02)/(1.25-0.75); 


%Plotting measurements points together with linear and theoretical fit
figure()
plot(Hrms_total, H13_total, 'ko', 'linewidth', 1)
hold on
plot(H_ref, fit, 'r-')
hold on 
plot(H_ref, theoretical_fit, 'b--')
legend('Measurements', 'Linear fit', 'Theoretical relation')
xlabel('H_{1/3} [m]')
ylabel('H_{rms} [m]')
title('Linear relation between H_{1/3} and H_{rms}')
xlim([0 1.6])

%% Relation between Hm and Hrms

%Use polyfit to determine the linear relation between Hrms and Hm
s = polyfit(Hrms_total, Hm_total, 1);
slope2 = s(1);
intercept2 = s(2);

%Linear fit using the slope and intercept found by polyfit
fit2 = slope2*H_ref + intercept2;

%Theoretical relation for mean and rms
Hm_theoretical = 0.89*H_ref; 

%Slope theoretical 
%0.25, 0.2225 
%1, 0.89

slope_theor2 = (0.89-0.225)/(1-0.25); 
%Plotting measurements points together with linear and theoretical fit 
figure()
plot(Hrms_total,Hm_total, 'ko', 'linewidth', 1)
hold on
plot(H_ref, fit2, 'r-')
hold on
plot(H_ref, Hm_theoretical, 'b--')
legend('Measurements', 'Linear fit', 'Theoretical relation')
xlabel('H_m [m]')
ylabel('H_{rms} [m]')
title('Linear relation between H_m and H_{rms}')

%% Save the different arrays containing the wave statistics for the full dataset

save('StatisticsEgmond', 'Hm_total', 'Hrms_total', 'H13_total')

