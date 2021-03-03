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

wlong = [-2.26 -3.27 4.68];

wcross = [8.86 8.55 1.16]; 

dzetady = [-7.33e-6 -1.76e-5 8.57e-7]; 

currents_data = load('vEgmond.txt'); 

position_sensors = [4765 4790 4814 4835 4860 4889]; %meters 

names = {'low tide', 'mid tide', 'high tide'}; 

for ii = 1:length(theta) 
     
waves = BJmodelEmma(Hrms0(ii),T0(ii),Zeta(ii),theta(ii),profile,hmin);


currents(:,ii) = longshoreCurrent(profile,dzetady(ii),wcross(ii),wlong(ii),waves.c,waves.theta,waves.Dr,waves.ht,waves.st,ka,nu,0)

%interpolated_v(:,ii) = interp1(waves.x,currents(:,ii),position_sensors);

figure() 
subplot(3,1,1)
plot(waves.x, waves.Hrms);
title(sprintf('Cross-shore evolution of Hrms for %s ',names{ii}));
grid on 
xlabel('x(m)') 
ylabel('Hrms(m)')
xlim([4000 max(waves.x)])
subplot(3,1,2)
plot(waves.x, currents(:,ii))
title('Modelled mean alongshore current and measured values') 
xlabel('x(m)') 
ylabel('v(m/s)')
xlim([4000 max(waves.x)])
grid on 
legend('Modelled along-shore current', 'Measured data'); 
hold on;
plot(position_sensors,currents_data(:,ii), '*') 
subplot(3,1,3) 
plot(bed_profile(:,1),bed_profile(:,2))
xlabel('x(m)') 
ylabel('z(m)') 
title('Bed level evolution') 
grid on 
     
end
 
%% Part 5.3: Complementary analysis 

%All analysis is for low tide only

for ii = 1:3 
    waves = BJmodelEmma(Hrms0(1),T0(1),Zeta(1),theta(1),profile,hmin);
    %curr(:,1) --> No wind forcing 
    %curr(:,2) --> No tidal forcing 
    %curr(:,3) --> No wind and no tidal forcing 
    curr(:,ii) = longshoreCurrent(profile,dzetady(1),wcross(1),wlong(1),waves.c,waves.theta,waves.Dr,waves.ht,waves.st,ka,nu,ii); 
    
end

%With waves arriving normal to the shore 

theta_normal = 90; 

waves_normal = BJmodelEmma(Hrms0(1),T0(1),Zeta(1),theta_normal,profile,hmin);


curr_normalwaves = longshoreCurrent(profile,dzetady(1),wcross(1),wlong(1),waves_normal.c,waves_normal.theta,waves_normal.Dr,waves_normal.ht,waves_normal.st,ka,nu,0)

figure()
subplot(3,1,1)  
plot(waves.x,currents(:,1)); 
title('Comparison of long-shore current with wind forcing and without forcing for low tide');
hold on; 
grid on; 
plot(waves.x,curr(:,1)); 
xlabel('x (m)'); 
ylabel('v (m/s)'); 
legend('With wind forcing', 'Without wind forcing'); 
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
plot(waves.x,currents(:,1)); 
title('Comparison of long-shore current with theta = -36 and theta = 90 for low tide'); 
hold on;
plot(waves_normal.x,curr_normalwaves); 
grid on; 
xlabel('x (m)'); 
ylabel('v (m/s)'); 
%legend('Theta = -36', 'Theta = 90'); 








 

 

 
 