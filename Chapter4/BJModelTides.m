%% Using the model with the Egmond tide data 

close all; 
clear all; 

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
     figure(1)
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
    figure(2) 
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
figure(1)
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














 

 
 
 
 

