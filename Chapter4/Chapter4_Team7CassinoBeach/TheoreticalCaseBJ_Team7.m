% Main script for the computation of the cross-shore wave evolution 
% for a simplified bar-trough case based on the Battjes and Janssen (1978) model 

%------------------------------------
%           Initialisation
%------------------------------------

clear all;
close all;

% Definition of cross-shore coordinates (m)
x = (1:1:500)';  

% Definition of zb, bed elevation relative to mean water level (m)
zb = -9 * ones(size(x));   
zb = zb + (x<=260).*(x-160)/20  + ((x<=300)&(x>260)).*(-1/40*(x-260) + 5) + (x>300).*(1/20*(x-300) + 4); 

% Definition of the array profile, input argument for BJmodel
profile = [x zb];

% Offshore wave conditions
Hrms0 = 1;      % Root mean square wave height (m)
theta0 = 0;     % Angle of incidence (degrees)
T0 = 10;        % Characteristic period (s)
Zeta = 0;       % Mean water level (m)

% Model parameter 
hmin = 0.2;     % Minimal water depth for computation
                % (we stop the computation when h<hmin)

%----------------------------------------------------------------------------------
%       Computation and visualisation of wave characteristics using a variable zeta
%----------------------------------------------------------------------------------

%% First case: Zeta = -1, 0 and 1 
Zeta = [-1 0 1]; %meters
theta0 = 0; %degrees 
Hrms0 = 1; %root mean square height 

%We run for all Zeta cases the BJ Model 
for i=1:length(Zeta);
    waves(i) = BJmodelEmma(Hrms0,T0,Zeta(i),theta0,profile,hmin);    
end 

% Case with Zeta = 0
figure(); 
%Modelled root mean square height
subplot(5,1,1); 
plot(waves(2).x, waves(2).Hrms);
ylabel('Hrms,0(m)');
xlim([0 500]);
ylim([0 2]); 
hold on;
title('Water level \zeta = 0m');
%Model set-up 
subplot(5,1,2); 
plot(waves(2).x, waves(2).eta);
ylabel('set-up (m)');
xlim([0 500]); 
ylim([-0.2 0.2]);
hold on;
%Dissipation of the breaking (Dbr) 
subplot(5,1,3);
plot(waves(2).x, waves(2).Dbr); 
ylabel('D_{br} (W/m^2)');
xlim([0 500]); ylim([0 400]);
hold on;
%Dissipation of the roller 
subplot(5,1,4);
plot(waves(2).x, waves(2).Dr); 
ylabel('D_{r} (W/m^2)');
xlim([0 500]);
ylim([0 400]); 
hold on;
%Bed profile 
subplot(5,1,5);
plot(waves(2).x, waves(2).z, 'k');
hold on;
plot(waves(2).x, Zeta(2)*ones(size(x)), '-.'); ylabel('zb (m)'); 
xlim([0 500]); ylim([-20 10]); xlabel('x(m)'); 
%% Figure with all cases of Zeta 


figure();
%Hrms for Zeta = -1 
subplot(5,1,1);
plot(waves(1).x,waves(1).Hrms) 
hold on; 
%Hrms for Zeta = 0  
plot(waves(2).x, waves(2).Hrms)
hold on; 
%Hrms for Zeta = 1 
plot(waves(3).x, waves(3).Hrms) 
legend('\zeta = -1', '\zeta =0', '\zeta = 1');
legend('location','northeast');
title('Water level \zeta = -1m, \zeta=0m and \zeta= 1m')
ylabel('Hrms,0 (m)')
xlim([0,500]); 
ylim([0 2]); 

%Plot of the setup 
subplot(5,1,2);
%Set up for Zeta = -1
plot(waves(1).x,waves(1).eta)
hold on; 
%Set up for Zeta = 0
plot(waves(2).x, waves(2).eta)
hold on;
%Set up for Zeta = 1
plot(waves(3).x, waves(3).eta)
xlim([0 500]); 
ylim([-0.2 0.2]);
ylabel('set-up (m)');

%Plot of Dbr 
subplot(5,1,3); 
%Dissipation of the breaking for Zeta = -1
plot(waves(1).x,waves(1).Dbr)
hold on; 
%Dissipation of the breaking for Zeta = 0
plot(waves(2).x, waves(2).Dbr)
hold on; 
%Dissipation of the breaking for Zeta = 1
plot(waves(3).x, waves(3).Dbr)
ylabel('D_{Br} (W/m^2)')
xlim([0,500]); 
ylim([0 400]);
%Plot of Dr 
subplot(5,1,4);
%Dissipation of the roller for Zeta = -1
plot(waves(1).x,waves(1).Dr)
hold on; 
%Dissipation of the roller for Zeta = 0
plot(waves(2).x, waves(2).Dr)
hold on; 
%Dissipation of the roller for Zeta = 1
plot(waves(3).x, waves(3).Dr)
xlim([0,500]); 
ylim([0 400]);
ylabel('D_r (W/m^2)')
%Plot of the bed level 
subplot(5,1,5); 
%Bed level with Zeta = -1
plot(waves(1).x,waves(1).z,'k')
hold on;
%Bed level with Zeta = 0
plot(waves(2).x, waves(2).z, 'k')
hold on; 
%Bed level with Zeta = 1
plot(waves(3).x, waves(3).z, 'k')
hold on; 
%Constant y line with Zeta = -1
plot(waves(1).x, Zeta(1)*ones(size(x)),'-.')
hold on; 
%Constant y line with Zeta = 0
plot(waves(2).x, Zeta(2)*ones(size(x)), '-.')
hold on; 
%Constant y line with Zeta = 1
plot(waves(3).x, Zeta(3)*ones(size(x)),'-.')
ylabel('zb (m)'); 
xlabel('x (m)');
xlim([0,500]); 
ylim([-20 10]);

%%
%----------------------------------------------------------------------------------
%       Computation and visualisation of wave characteristics using a
%       variable Hrms
%----------------------------------------------------------------------------------


Zeta = 0 %meters
Hrms0 = [0.5, 2]; %meters
theta0=0;
Title = ["Hrms0 = 0.5m", "Hrms0 = 2m"];

%We run the cases for the two different values of Hrms0 
for i=1:length(Hrms0)
    waves(i) = BJmodelEmma(Hrms0(i),T0,Zeta,theta0,profile,hmin);
end 


%Figure of wave transformation with different values of Hrms0 
figure();
%Plot of Hrms 
subplot(5,1,1); 
% Plot of modelled Hrms with Hrms0 = 0.5 
plot(waves(1).x,waves(1).Hrms)
hold on; 
% Plot of modelled Hrms with Hrms0 = 2
plot(waves(2).x, waves(2).Hrms)
title('Hrms = 0.5m and 2 m, \zeta =0m')
ylabel('Hrms,0 (m)'); 
legend('Hrms=0.5', 'Hrms=2')
xlim([0,500]); 
ylim([0 3])
%Plot of the bed level from BJModel 
subplot(5,1,5);
%Plot of bed level with Hrms0 = 0.5 
plot(waves(1).x,waves(1).z,'k')
hold on; 
%Plot of bed level with Hrms0 = 2
plot(waves(2).x, waves(2).z, 'k')
hold on; 
%Constant line with Zeta = 0 
plot(waves(1).x,Zeta*ones(size(x)),'-.')
ylabel('zb (m)'); 
xlabel('x (m)');
xlim([0,500]);
%Dissipation of the breaking 
subplot(5,1,3);
%Dissipation of the breaking for Hrms0 = 0.5 
plot(waves(1).x,waves(1).Dbr)
hold on; 
%Dissipation of the breaking for Hrms0 = 2
plot(waves(2).x, waves(2).Dbr)
ylabel('D_{Br} (W/m^2)')
xlim([0,500]);
%Dissipation of the roller
subplot(5,1,4); 
%Dissipation of the roller for Hrms0 = 0.5 
plot(waves(1).x,waves(1).Dr)
hold on; 
%Dissipation of the roller for Hrms0 = 2
plot(waves(2).x, waves(2).Dr)
xlim([0,500]);
ylabel('D_r (W/m^2)')
%Plot of the setup 
subplot(5,1,2); 
%Set up with Hrms = 0.5
plot(waves(1).x,waves(1).eta)
hold on; 
%Set up with Hrms = 2
plot(waves(2).x, waves(2).eta)
xlim([0,500]); 
ylim([-0.2 0.4])
ylabel('set-up (m)')


%----------------------------------------------------------------------------------
%       Computation and visualisation of wave characteristics using a
%       variable theta
%----------------------------------------------------------------------------------
%% Now we run the case with different values for theta (the incidence angle) 
Zeta = 0 %meters
Hrms0 = 1 %meters
theta0 = [22.5, 45]; %degrees
Title = ["theta0 = 22.5 degrees", "theta0 = 45 degrees"];

%We run the BJModel for both cases of theta 
for i=1:length(theta0);
    waves(i) = BJmodelEmma(Hrms0,T0,Zeta,theta0(i),profile,hmin);
end


%Modelled Hrms
figure();
subplot(5,1,1);
%Modelled Hrms with theta = 22.5 
plot(waves(1).x,waves(1).Hrms)
hold on;
%Modelled Hrms with theta = 45
plot(waves(2).x, waves(2).Hrms)
title('\theta = 22.5^{\circ} and \theta=45^{\circ}') ; 
legend('\theta = 22.5^{\circ}', '\theta = 45^{\circ}');
ylabel('Hrms (m)')
xlim([0,500]); 
ylim([0 2]);
%Plot of set up 
subplot(5,1,2); 
%Set up with theta = 22.5 
plot(waves(1).x,waves(1).eta)
hold on; 
%Set up with theta = 45 
plot(waves(2).x, waves(2).eta)
xlim([0,500]); 
ylim([-0.2 0.2])
ylabel('set-up (m)')
%Plot of dissipation of the breaking 
subplot(5,1,3); 
%Dissipation of the breaking with theta = 22.5 
plot(waves(1).x,waves(1).Dbr)
hold on; 
%Dissipation of the breaking with theta = 45 
plot(waves(2).x, waves(2).Dbr)
ylabel('D_{Br} (W/m^2)')
xlim([0,500]);
ylim([0 400]);
%Plot of the dissipation of the roller 
subplot(5,1,4); 
%Dissipation of the roller for theta = 22.5
plot(waves(1).x,waves(1).Dr)
hold on; 
%Dissipation of the roller for theta = 45
plot(waves(2).x, waves(2).Dr)
xlim([0,500]);
ylim([0 400]);
ylabel('D_r (W/m^2)')
%Plot of the bed level 
subplot(5,1,5); 
%Bed level for theta = 22.5 
plot(waves(1).x,waves(1).z,'k')
hold on; 
%Bed level for theta = 45
plot(waves(2).x, waves(2).z, 'k')
hold on; 
%Constant y line with Zeta = 0 
plot(waves(1).x,Zeta*ones(size(x)),'-.')
hold on; 
ylabel('zb (m)'); 
xlabel('x (m)')
xlim([0,500]);


%% Using the model with the Egmond tide data 

% Values for theta, T and H1/3 from Table 9.1 

%Low tide, mid tide and high tide 

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

%We load the statistics for the Egmond data that we calculate for Chapter 1
load('StatisticsEgmond.mat') 

%first column: positions (x), second column: water level (z) 
profile = [bed_profile(:,1) bed_profile(:,2)]; 


% Model parameter 
hmin = 0.2;     % Minimal water depth for computation
                % (we stop the computation when h<hmin)
                
%This names are for the title for each one of the tide plots 
names = {'low tide', 'mid tide', 'high tide'};

%This are the positions of each one of the sensors for Egmond 
position_sensors = [4478 4765 4790 4814 4835]; %meters 
                
 
%We run the for loop to calculate the BJ Model for each tide 
 for ii  = 1:length(theta)
     %We calculate the BJ Model with respect of each tide Hrms, period,
     %angle of incidence
     %All tides use the same profile and same h min 
     waves = BJmodelEmma(Hrms0(ii),T0(ii),Zeta(ii),theta(ii),profile,hmin);
   
     %We create the vertical subplots showing the modelled Hrms and
     %comparing it to the observed Hrms calculated on Chapter 1 for each
     %tide 
     figure(5)
     %We iterate over the different vertical subplots to plot each tide 
     subplot(4,1,ii)
     %We plot the modelled Hrms 
     plot(waves.x,waves.Hrms)
     %Title of the plot 
     title(sprintf('Comparison of modelled Hrms and observed Hrms for %s ',names{ii})); 
     hold on; 
     %We create the colors to do a plot of the Hrms observations (Chapter
     %1) 
     cmap = winter(length(position_sensors)); % Make colors 
     for jj= 1:length(position_sensors)
     scatter(position_sensors(jj), Hrms_total(jj,ii), 10, cmap(jj,2), 'filled')
     leg = {'Hrms BJModel', 'Hrms:P1','Hrms:P2','Hrms:P3','Hrms:P4','Hrms:P5'};  
     %We define the labels and the range for the x and y axis 
     ylabel('Hrms (m)')
     xlabel('x (m)') 
     xlim([4000 5000])
     ylim([0.3 2]) 
     grid on 
     end
     legend(leg)
 
     
    %Interpolation of the modelled Hrms into the positions of the sensors 
    interpolation_data{ii} = interp1(waves.x,waves.Hrms,position_sensors);
    
    %Figure to check that the interpolation of the modelled Hrms is done
    %correctly 
    figure(6) 
    subplot(3,1,ii) 
    %Plot of the interpolated Hrms 
    plot(position_sensors,interpolation_data{ii}, '*')
    title(sprintf('Interpolated Hrms model values vs observed Hrms values for %s', names{ii}));  
    hold on 
    %Plot of the observed Hrms (Chapter 1) 
    plot(position_sensors,Hrms_total(:,ii),'*')
    legend('Hrms model', 'Hrms observations') 
    xlabel('x(m)') 
    ylabel('Hrms(m)') 
    ylim([0.2 1.7]) 
    grid on 
    
    %Calculation of the root mean square error for each tide 
    %Xobs --> Hrms observations (Chapter 1) 
    %Xmodel --> Interpolated Hrms 
    rmse(ii) = rootmeansquare_error(Hrms_total(:,ii),interpolation_data{ii}); 
     
 end

 
%Last subplot of Figure 1 with the bed evolution and determining regions
%for each one of the tides 
figure(5)
subplot(4,1,4)
hold on;
%Plot of the bed profile 
plot(bed_profile(:,1),bed_profile(:,2))
title('Bed level evolution') 
hold on;
%Plot of y lines with water level for each tide
%Low tide 
yline(Zeta(1), '--','color', 'blue') 
%Mid tide 
yline(Zeta(2),'--','color', 'magenta') 
%High tide 
yline(Zeta(3),'--','color', 'black') 
legend('Bed level', 'low tide', 'mid tide', 'high tide') 
grid on 
xlim([4000 max(bed_profile(:,1))])
ylabel('z (m)') 
xlabel('x (m)') 
ylim([-10 3])


% Calculating the root mean square error for the whole series 

Obs_together = [transpose(Hrms_total(:,1)) transpose(Hrms_total(:,2)) transpose(Hrms_total(:,3))]; 

Model_together = [interpolation_data{1} interpolation_data{2} interpolation_data{3}]; 

rmse_allseries = rootmeansquare_error(Obs_together,Model_together);




