%% Chapter 6 

clear all; 
close all; 



%% 6.1.1 Computation from see surface elevation time series 

%We first declare the constants from Table 9.1 that we will use when
%computing waves 

%low tide, mid tide and high tide 

%Angle of wave incidence 
theta = [-36 39 36]; %degrees 

%Period 1/3 
T0 = [7.58 6.69 5.54]; %seconds 

%Significant wave height 
H13 = [1.70 2.25 1.69]; %meters 

%Calculated root mean square wave height (as calculated H1/3 =
%sqrt(2)*Hrms)
Hrms0 = H13./sqrt(2); %meters 

%Water level (m NAP) 
Zeta = [-0.45 0.09 0.91] %meters 


%Bed profile for Egmond data 
bed_profile = load('prof1018.txt'); 

%first column: positions (x), second column: water level (z) 
profile = [bed_profile(:,1) bed_profile(:,2)]; 


% Model parameter 
hmin = 0.2;     % Minimal water depth for computation
                % (we stop the computation when h<hmin)
                

%We load the mean water depth (it will be h) 
                
h = load('MeanWaterDepth.txt'); 


%Position of the sensors 
position_sensors = [4478 4765 4790 4814 4835]; %meters 


array_names = ["lowTide.txt", "midTide.txt", "highTide.txt"];

%We load the statistics from Chapter 1 

load('StatisticsEgmond.mat');

%We initialize the matrices that will contain our skewness and assymetry as
%well as the Ursell number 
skewness = zeros(5,3); 
assymetry = zeros(5,3); 
Ur = zeros(5,3); 

%First we run for all position of the sensors 

for jj = 1:length(position_sensors)
    
    %We run for each one of the tides 
    
for ii = 1:length(array_names) 
    
    %We iterate the name of the file (low tide, mid tide or high tide) 
    name = array_names(ii);
    wave_data = load(name);
    
    wavenumber = wavenumber_Guo(T0(ii), h(jj,ii));  
    %We calculate the skewness and assymetry for each one of the tides for
    %all sensors 
    [skewness(jj,ii) assymetry(jj,ii)] = skewness_asymmetry(wave_data(:,jj));
    %We calculate the ursell number for each one of the tides for all
    %sensors 
    Ur_obs(jj,ii) = ursell_number(wavenumber,h(jj,ii),Hrms_total(jj,ii)); 
    
     
end  

   
end

%We save the skewness, assymetry and ursell number for all tides for all
%sensors 
save('Skewness_assymetry_Ur.mat', 'skewness','assymetry','Ur_obs');  


%% 6.1.2 Comparison with the empirical fit defined by Ruessink et al (2012)

% Vector of Ur varying from 0.01 to 100

Ur = [0.01:0.01:100]; 

%calculation of skewness and assymetry as Ruessink for our Ur vector 

[skewness_R assymetry_R] = skewness_assymetryRuessink(Ur); 

%Plot of skweness and assymetry of Ruessink vs previous skewness and
%assymetry for all tides all sensors 
figure(1)
subplot(2,1,1) 
%Plot skewness Ruessink 
semilogx(Ur,skewness_R) 
hold on;
%Plot skewness low tide 
semilogx(Ur_obs(:,1),skewness(:,1),'*'); 
%Plot skewness mid tide 
semilogx(Ur_obs(:,2),skewness(:,2),'*'); 
%Plot skewness high tide 
semilogx(Ur_obs(:,3),skewness(:,3),'*'); 
title('Skewness as a function of the Ursell number')
xlabel('Ursell number (log scale)')
ylabel('Skewness') 
grid on 
legend('Skewness Ruessink', 'Skewness low tide', 'Skewness mid tide', 'Skewness high tide'); 
%Plot assymetry Ruessink and tides 
subplot(2,1,2)
%Plot assymetry Ruessink 
semilogx(Ur,assymetry_R) 
hold on
%Plot assymetry low tide 
semilogx(Ur_obs(:,1),assymetry(:,1),'*'); 
%Plot assymetry mid tide 
semilogx(Ur_obs(:,2),assymetry(:,2),'*'); 
%Plot assymetry high tide 
semilogx(Ur_obs(:,3),assymetry(:,3),'*'); 
title('Assymetry as a function of the Ursell number')
xlabel('Ursell number (log scale)')
ylabel('Assymetry') 
legend('Assymetry Ruessink', 'Assymetry low tide', 'Assymetry mid tide', 'Assymetry high tide'); 
grid on 

%% 6.2 Modelling Sk and As 

title_it = {'low tide', 'mid tide', 'high tide'}; 

%Colors for each one of the sensors (P1 P3 P4 P5 P6) 
c = {'r', 'g', 'b', 'm', 'c'}; 

%We iterate for each one of the tides 
for ii = 1: length(theta) 
%We get eta from the BJModel 
waves = BJmodelEmma(Hrms0(ii),T0(ii),Zeta(ii),theta(ii),profile,hmin);


%We now calculate the ursell number 
ur_waves = ursell_number(waves.k, waves.ht, waves.Hrms);

%Now we calculate As and Sk with the empirical relations 

[sk_waves as_waves] = skewness_assymetryRuessink(ur_waves);


%We plot the skewness calculated with waves for each tide 
figure() 
subplot(3,1,1)
plot(waves.x, sk_waves); 
title(sprintf('Cross-shore evolution of the skewness predicted by the BJModel for %s', title_it{ii})); 
hold on 
for jj = 1:length(position_sensors)
scatter(position_sensors(jj),skewness(jj,ii),10, c{jj}, 'filled')
leg = {'Sk:Ruessink', 'Sk:P1','Sk:P3','Sk:P4','Sk:P5','Sk:P6'};  
end
legend(leg)
xlabel('x(m)') 
ylabel('Skewness') 
grid on 
xlim([4000 max(waves.x)]); 

%Plot of assymetry for each tide 
subplot(3,1,2) 
plot(waves.x, as_waves)
title(sprintf('Cross-shore evolution of the assymetry predicted by the BJModel for %s', title_it{ii})); 
hold on 
for jj = 1:length(position_sensors)
scatter(position_sensors(jj),assymetry(jj,ii),10, c{jj}, 'filled')
leg = {'As:Ruessink', 'As:P1','As:P3','As:P4','As:P5','As:P6'};  
end
legend(leg)
grid on 
xlabel('x(m)') 
xlim([4000 max(waves.x)]); 
ylabel('Assymetry') 
%Bed level evolution 
subplot(3,1,3)
plot(bed_profile(:,1),bed_profile(:,2))
xlabel('x(m)') 
ylabel('z(m)') 
xlim([4000 max(bed_profile(:,1))]); 
ylim([min(bed_profile(:,2)) max(bed_profile(:,2))]); 
title('Bed level evolution') 
grid on 



end



%% 6.3 Cross-shore evolution of the orbital velocity 

%Phi for all cases 
phi = [-pi/2 -pi/4 0];

%R for all cases 
r = [0 0.3 0.6];

%Period 
T = 10 %seconds 
%Velocity amplitude 
Uw = 1 %m/s 

%Color for colorbar 
C = {'k', 'm', 'b'};

%We run for all cases 
for jj = 1:length(phi)
%First case : phi = -pi/2 and r =0, 0.3, 0.6 
    
[u{jj} t{jj}] = waveshape(r(jj),phi(1),Uw,T); 

figure(5)
subplot(3,1,1) 
hold on;
plot(t{jj},u{jj})
title('Orbital velocity with \phi = - \pi /2'); 
xlabel('t (s)')
ylabel(' u (m/s)')
grid on; 

%Second case: phi = 0 and r = 0 0.3, 0.6 

[u2{jj} t2{jj}] = waveshape(r(jj),phi(3),Uw,T);

subplot(3,1,2)
hold on; 
plot(t2{jj},u2{jj})
title('Orbital velocity with \phi = 0'); 
xlabel('t (s)')
ylabel(' u (m/s)')
grid on; 

%Third case: r = 0.6 and phi = 0, -pi/2, -pi/4 

[u3{jj} t3{jj}] = waveshape(r(3),phi(jj),Uw,T);

subplot(3,1,3)
hold on; 
plot(t3{jj},u3{jj})%,'color', C{ii})
title('Orbital velocity with r = 0.6'); 
xlabel('t (s)')
ylabel(' u (m/s)')
grid on; 


end

%We add legends to the plot 
figure(5)
subplot(3,1,1)
legend('r = 0', 'r = 0.3','r = 0.6'); 
subplot(3,1,2)
legend('r = 0', 'r = 0.3','r = 0.6'); 
subplot(3,1,3) 
legend(' \phi = - \pi /2 ', ' \phi = - \pi /4 ', ' \phi = 0 '); 

%% Low tide case 

%We calculate waves for the low tide     
waves = BJmodelEmma(Hrms0(1),T0(1),Zeta(1),theta(1),profile,hmin);

%We calculate velocity amplitude for the low tide     
Uw_lowtide = velocity_amplitude(waves.ht,waves.Hrms,T0(1));

%Plot of the evolution of the velocity amplitude for the low tide 
figure()
subplot(2,1,1)
plot(waves.x, Uw_lowtide)
hold on
grid on 
title('Cross-shore evolution of the velocity amplitude for low tide')
xlabel('z(m)')
xlim([4000 5000]) 
ylim([0 1]) 
ylabel('Velocity amplitude (m)')
subplot(2,1,2) 
plot(bed_profile(:,1),bed_profile(:,2), 'k')
xlabel('x(m)') 
ylabel('z(m)') 
xlim([4000 max(bed_profile(:,1))]); 
ylim([min(bed_profile(:,2)) max(bed_profile(:,2))]); 
title('Bed level evolution') 
grid on 

%% Different positions 

x = [1000 4400 4500 4700 4920]; %positions (meters) 


[pos, ~] = find(waves.x == x); % find indices of x-positions

%Now we calculate As and Sk with the empirical relations 

%Ursell number calculations with output of waves  
ur_waves_lowtide = ursell_number(waves.k, waves.ht, waves.Hrms);

%Skewness and assymetry with output of waves 

[sk_lowtide as_lowtide] = skewness_assymetryRuessink(ur_waves_lowtide);



%Computation of r and phi for the low tide and the new positions

for xx = 1:length(pos) 

%Calculation r 
r_lowtide(xx) = computation_r(sk_lowtide(pos(xx,1)),as_lowtide(pos(xx,1))); 

%Calculation phi 
phi_lowtide(xx) = computation_phi(sk_lowtide(pos(xx,1)),as_lowtide(pos(xx,1))); 

%Velocity amplitude 
Uw_newpositions(xx) = Uw_lowtide(pos(xx,1)); 

%Orbital velocity 
[u_new(:,xx) t_new(:,xx)] = waveshape(r_lowtide(xx),phi_lowtide(xx),Uw_newpositions(xx),T0(1)); 


%Time series of orbital velocity for each one of the positions 
figure(10)
subplot(5,1,xx)
plot(t_new(:,xx),u_new(:,xx))
xlabel('t(s)') 
ylabel(' u(t) (m/s)')
title(sprintf('Orbital velocity for position x = %d m',x(xx)))
grid on 
ylim([-1.5 1.5]) 


end








