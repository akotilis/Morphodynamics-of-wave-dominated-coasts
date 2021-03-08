%% Chapter 6 


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


array_names = ["lowTide.txt", "midTide.txt", "highTide.txt"];

skewness = zeros(5,3); 
assymetry = zeros(5,3); 
Ur = zeros(5,3); 

for jj = 1:length(wave_data(1,:))
    
for ii = 1:length(names_datasets) 
    
    name = array_names(ii);
    wave_data = load(name);
    
    %wavenumber = wavenumber_Guo(T0(ii), h) 
    
    [skewness(jj,ii) assymetry(jj,ii)] = skewness_asymmetry(wave_data(:,jj));
    
    %Ur(ii) = ursell_number(k,h,Hrms)
    
    %save(sprintf('Skewness_%s', name)
    
end  
   
end


%% 6.1.2 

% Vector of Ur varying from 0.01 to 100

Ur = [0.01:0.01:100]; 

[skewness_R assymetry_R] = skewness_assymetryRuessink(Ur); 

figure()
subplot(2,1,1) 
semilogx(Ur,skewness_R) 
title('Skewness as a function of the Ursell number')
xlabel('Ursell number (log scale)')
ylabel('Skewness') 
grid on 

subplot(2,1,2) 
semilogx(Ur,assymetry_R) 
title('Assymetry as a function of the Ursell number')
xlabel('Ursell number (log scale)')
ylabel('Assymetry') 
grid on 

%% 6.2 Modelling Sk and As 

%We first get the eta from the BJModel 

%Bed profile for Egmond data 
bed_profile = load('prof1018.txt'); 

%first column: positions (x), second column: water level (z) 
profile = [bed_profile(:,1) bed_profile(:,2)]; 


% Model parameter 
hmin = 0.2;     % Minimal water depth for computation
                % (we stop the computation when h<hmin)

for ii = 1: length(theta) 
    %We get eta from the BJModel 
waves = BJmodelEmma(Hrms0(ii),T0(ii),Zeta(ii),theta(ii),profile,hmin);


%We now calculate the ursell number 


ur = ursell_number(waves.k, waves.eta, waves.Hrms);

figure() 
subplot(3,1,1) 

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







