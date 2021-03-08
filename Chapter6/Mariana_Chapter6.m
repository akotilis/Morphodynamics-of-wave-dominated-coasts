%% Chapter 6 

clear all; 
close all; 



%% 6.1.1 

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


%We first get the eta from the BJModel 

%Bed profile for Egmond data 
bed_profile = load('prof1018.txt'); 

%first column: positions (x), second column: water level (z) 
profile = [bed_profile(:,1) bed_profile(:,2)]; 


% Model parameter 
hmin = 0.2;     % Minimal water depth for computation
                % (we stop the computation when h<hmin)
                
h = load('MeanWaterDepth.txt'); 

position_sensors = [4478 4765 4790 4814 4835]; %meters 


array_names = ["lowTide.txt", "midTide.txt", "highTide.txt"];

skewness = zeros(5,3); 
assymetry = zeros(5,3); 
Ur = zeros(5,3); 

for jj = 1:length(wave_data(1,:))
    
for ii = 1:length(array_names) 
    
    name = array_names(ii);
    wave_data = load(name);
    
    wavenumber = wavenumber_Guo(T0(ii), h(jj,ii));  
    
    [skewness(jj,ii) assymetry(jj,ii)] = skewness_asymmetry(wave_data(:,jj));
    
    Ur_obs(jj,ii) = ursell_number(wavenumber,h(jj,ii),Hrms0(ii)); 
    
     
end  

   
end

save('Skewness_assymetry_Ur.mat', 'skewness','assymetry','Ur_obs');  


%% 6.1.2 

% Vector of Ur varying from 0.01 to 100

Ur = [0.01:0.01:100]; 

[skewness_R assymetry_R] = skewness_assymetryRuessink(Ur); 

figure()
subplot(2,1,1) 
semilogx(Ur,skewness_R) 
hold on;
semilogx(Ur_obs(:,1),skewness(:,1),'*'); 
semilogx(Ur_obs(:,2),skewness(:,2),'*'); 
semilogx(Ur_obs(:,3),skewness(:,3),'*'); 
title('Skewness as a function of the Ursell number')
xlabel('Ursell number (log scale)')
ylabel('Skewness') 
grid on 
legend('Skewness Ruessink', 'Skewness low tide', 'Skewness mid tide', 'Skewness high tide'); 

subplot(2,1,2) 
semilogx(Ur,assymetry_R) 
hold on
semilogx(Ur_obs(:,1),assymetry(:,1),'*'); 
semilogx(Ur_obs(:,2),assymetry(:,2),'*'); 
semilogx(Ur_obs(:,3),assymetry(:,3),'*'); 
title('Assymetry as a function of the Ursell number')
xlabel('Ursell number (log scale)')
ylabel('Assymetry') 
legend('Assymetry Ruessink', 'Assymetry low tide', 'Assymetry mid tide', 'Assymetry high tide'); 
grid on 

%% 6.2 Modelling Sk and As 

title_it = {'low tide', 'mid tide', 'high tide'}; 


for ii = 1: length(theta) 
%We get eta from the BJModel 
waves = BJmodelEmma(Hrms0(ii),T0(ii),Zeta(ii),theta(ii),profile,hmin);


%We now calculate the ursell number 


ur_waves = ursell_number(waves.k, waves.eta, waves.Hrms);

%Now we calculate As and Sk with the empirical relations 

[sk_waves as_waves] = skewness_assymetryRuessink(ur_waves);


figure() 
subplot(3,1,1)
plot(waves.x, sk_waves(ii)); 
title(sprintf('Cross-shore evolution of the skewness predicted by the BJModel for %s', title_it{ii})); 
hold on 
semilogx(Ur_obs(:,ii),skewness(:,ii),'*') 
legend('Skewness with Ruessink', 'Skewness observations')
xlabel('Ursell number (log scale)')
ylabel('Skewness') 
grid on 

subplot(3,1,2) 
plot(ur_waves, as_waves(ii))
title(sprintf('Cross-shore evolution of the assymetry predicted by the BJModel for %s', title_it{ii})); 
hold on 
plot(Ur_obs(:,ii),assymetry(:,ii),'*'); 
legend('Assymetry with Ruessink', 'Assymetry observations')
grid on 
xlabel('Ursell number (log scale)')
ylabel('Assymetry') 

subplot(3,1,3)
plot(bed_profile(:,1),bed_profile(:,2))
xlabel('x(m)') 
ylabel('z(m)') 
xlim([0 max(bed_profile(:,1))]); 
ylim([min(bed_profile(:,2)) max(bed_profile(:,2))]); 
title('Bed level evolution') 
grid on 



end



%% 6.3 Cross-shore evolution of the orbital velocity 

phi = [-pi/2 -pi/4 0];

r = [0 0.3 0.6];

T = 10 %seconds 
Uw = 1 %m/s 

C = {'k', 'm', 'b'};

for jj = 1:length(phi)
%First case : phi = -pi/2 and r =0, 0.3, 0.6 
    
[u{jj} t{jj}] = waveshape(r(jj),phi(1),Uw,T); 

figure(1)
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

%Plots 
figure(1)
subplot(3,1,1)
legend('r = 0', 'r = 0.3','r = 0.6'); 
subplot(3,1,2)
legend('r = 0', 'r = 0.3','r = 0.6'); 
subplot(3,1,3) 
legend(' \phi = - \pi /2 ', ' \phi = - \pi /4 ', ' \phi = 0 '); 



%Now we plot for the low tide case 


for xx = 1:length(h(:,1))
    
Uw(xx) = velocity_amplitude(h(xx,1),Hrms0(1),T0(1)) 

end

figure()
plot(position_sensors,Uw, '*')
grid on 
title('Cross-shore evolution of the velocity amplitude for low tide')
xlabel('z(m)')
ylabel('Velocity amplitude (m)')


%Now we calculate for different positions 



x = [1000 4400 4500 4700 4920]; %positions (meters) 

% It takes skewness, assymetry computation_r(S,A)

% computation_phi(S,A)









