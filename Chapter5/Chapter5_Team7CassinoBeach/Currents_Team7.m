%% Chapter 5 

%% Section 5.2: Data/Model Comparison for Egmond 
close all; 
clear all; 
%We first put the initial conditions for each tide to run the BJ Model 

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


%Constants


ka = 0.022 % Apparent bed roughness [meters]

nu =  0.5  %large-scale mixing coefficient [mˆ2/s]

%Constants from Table 9.1 


%Along-shore wind speed 
wlong = [-2.26 -3.27 4.68]; %m/s 

%Cross-shore wind speed 
wcross = [8.86 8.55 1.16]; %m/s 


dzetady = [-7.33e-6 -1.76e-5 8.57e-7]; 

%We load the observations for the current for Egmond 
currents_data = load('vEgmond.txt'); 
%'P1','P3','P4','P5','P6'
position_pastsensors = [4478 4765 4790 4814 4835]; %meters 

%'P3', 'P4', 'P5', 'P7' 'P8'
%This are the position of the sensors for the current data 
position_sensors = [4765 4790 4814 4835 4860 4889]; %meters 

%We load the statistics from Chapter 1 
load('StatisticsEgmond.mat'); 

names = {'low tide', 'mid tide', 'high tide'}; 

c = {'r', 'g', 'b', 'm', 'c'}; 

c2 = {'g', 'b', 'm', 'c', 'k', 'y'};

%We iterate for each one of the tides 
for ii = 1:length(theta) 

%We calculate the wave properties over the profile for each tide 
waves = BJmodelEmma(Hrms0(ii),T0(ii),Zeta(ii),theta(ii),profile,hmin);


%We calculate the current for each one of the tides 
currents(:,ii) = longshoreCurrent(profile,dzetady(ii),wcross(ii),wlong(ii),waves.c,waves.theta,waves.Dr,waves.ht,waves.st,ka,nu,0);

%Interpolation of the current profile to the position of the sensors 
interpolated_v{ii} = interp1(waves.x,currents(:,ii),position_sensors);

%Calculation of root mean square error for each tide 
rmse(ii) = rootmeansquare_error(currents_data(:,ii),interpolated_v{ii});
sprintf('The root mean square error for %s is %d', names{ii}, rmse(ii))

figure() 
%Plot of modelled Hrms 
subplot(3,1,1)
plot(waves.x, waves.Hrms);
title(sprintf('Cross-shore evolution of Hrms for %s ',names{ii}));
grid on 
ylim([min(waves.Hrms) max(waves.Hrms)]); 
hold on; 
for jj= 1:length(position_pastsensors)
scatter(position_pastsensors(jj), Hrms_total(jj,ii), 10, c{jj}, 'filled')
leg = {'Hrms BJModel', 'Hrms:P1','Hrms:P3','Hrms:P4','Hrms:P5','Hrms:P6'};  
%We define the labels and the range for the x and y axis 
end
legend(leg)
xlabel('x(m)') 
ylabel('Hrms(m)')
xlim([4000 max(waves.x)])
%Plot of the cross-shore distribution of the longshore current
subplot(3,1,2)
plot(waves.x, currents(:,ii))
hold on;
title('Modelled mean alongshore current and measured values') 
xlabel('x(m)') 
ylabel('v(m/s)')
xlim([4000 max(waves.x)])
%ylim([min(currents(:,ii)) max(currents(:,ii))]); 
grid on 
for xx= 1:length(position_sensors)
scatter(position_sensors(xx), currents_data(xx,ii), 10, c2{xx}, 'filled')
lege = {'Modelled along-shore current', 'Hrms:P3','Hrms:P4','Hrms:P5','Hrms:P7','Hrms:P8'};  
%We define the labels and the range for the x and y axis 
end
legend(lege)
%Plot of the bed profile 
subplot(3,1,3) 
plot(bed_profile(:,1),bed_profile(:,2))
xlabel('x(m)') 
ylabel('z(m)') 
xlim([0 max(bed_profile(:,1))]); 
ylim([min(bed_profile(:,2)) max(bed_profile(:,2))]); 
title('Bed level evolution') 
grid on 
     
end


 
%% Part 5.3: Complementary analysis 

%All analysis is for low tide only

%We added an if statement in the function longshoreCurrent.m to run for 
%the different cases that the question asks for, defining an input type,
%and depending on the number is the forcings it takes into account for the
%F calculation 

for ii = 1:3 
    waves = BJmodelEmma(Hrms0(1),T0(1),Zeta(1),theta(1),profile,hmin);
    %curr(:,1) --> No wind forcing 
    %curr(:,2) --> No tidal forcing 
    %curr(:,3) --> No wind and no tidal forcing 
    curr(:,ii) = longshoreCurrent(profile,dzetady(1),wcross(1),wlong(1),waves.c,waves.theta,waves.Dr,waves.ht,waves.st,ka,nu,ii); 
    
end

%With waves arriving normal to the shore 

theta_normal = 0; 

waves_normal = BJmodelEmma(Hrms0(1),T0(1),Zeta(1),theta_normal,profile,hmin);


curr_normalwaves = longshoreCurrent(profile,dzetady(1),wcross(1),wlong(1),waves_normal.c,waves_normal.theta,waves_normal.Dr,waves_normal.ht,waves_normal.st,ka,nu,0);

figure()
%Plot of long-shore current with wind forcing and without wind forcing 
subplot(3,1,1)  
plot(waves.x,currents(:,1)); 
title('Comparison of long-shore current with wind forcing and without forcing for low tide');
hold on; 
grid on; 
plot(waves.x,curr(:,1)); 
xlabel('x (m)'); 
ylabel('v (m/s)'); 
legend('With wind forcing', 'Without wind forcing'); 
%Plot of comparison of long-shore current without wind, and without wind and tidal forcing
%
subplot(3,1,2)
plot(waves.x,curr(:,1)); 
title('Comparison of long-shore current without wind and tidal forcing and without tidal forcing for low tide'); 
hold on; 
grid on; 
plot(waves.x,curr(:,3)); 
xlabel('x (m)'); 
ylabel('v (m/s)'); 
legend('Without wind forcing', 'Without wind and without tidal forcing') 
subplot(3,1,3) 
%Plot of long-shore current with theta = -36 and waves normal to the sea
%shore 
plot(waves.x,currents(:,1)); 
title('Comparison of long-shore current with theta = -36 and theta = 0 for low tide'); 
hold on;
plot(waves_normal.x,curr_normalwaves); 
grid on; 
xlabel('x (m)'); 
ylabel('v (m/s)'); 
legend('Theta = -36', 'Theta = 0'); 








 

 

 
 