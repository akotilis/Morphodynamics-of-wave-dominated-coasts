%% 8.2  Sediment transport by waves only 
close all; 
clear all; 


%Calculation of orbital velocity for mid tide 

%Angle of wave incidence 
theta = 39; %degrees 

%Period 1/3 
T0 = 6.69; %seconds 

%Significant wave height 
H13 = 2.25; %meters 

%Calculated root mean square wave height (as calculated H1/3 =
%sqrt(2)*Hrms)
Hrms0 = H13./sqrt(2); %meters 

%Water level (m NAP) 
Zeta = 0.09; %meters 


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


%array_names = "midTide.txt"
   
waves = BJmodelEmma(Hrms0,T0,Zeta,theta,profile,hmin);

%We calculate velocity amplitude for the low tide     
Uw = velocity_amplitude(waves.ht,waves.Hrms,T0);
 
ur_waves = ursell_number(waves.k, waves.ht, waves.Hrms);

%Skewness and assymetry with output of waves 

[sk as] = skewness_assymetryRuessink(ur_waves);

N_last = find(~isnan(sk),1,'last');

%Computation of r and phi for mid tide 

for xx = 1:N_last
    
%Calculation r 
r(xx) = computation_r(sk(xx,1),as(xx,1)); 

%Calculation phi 
phi(xx) = computation_phi(sk(xx,1),as(xx,1)); 

%Orbital velocity 
%[u(:,xx) t(:,xx)] = waveshape(r(xx),phi(xx),Uw(xx),T0); 
end 

save('BattjesJansenmidtideChapter6.mat', 'waves', 'Uw', 'ur_waves', 'sk', 'as', 'r', 'phi', 'N_last')
