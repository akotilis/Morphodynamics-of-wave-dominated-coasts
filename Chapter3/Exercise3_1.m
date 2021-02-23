%% Chapter 3 : Linear wave theory 

clear all; 
close all; 

periods = [6 9 12]; %seconds 
water_depth = [0:0.5:140]; %meters

%Constants 

g  = 9.81; %Gravity, m/sˆ2 
omega = 2*pi./periods;

%We run a for loop for each one of the periods 
for ii = 1:length(periods)
    %We calculate the wavenumber with the function we made, taking our
    %water depth vector and the period we are iterating on 
wavenumber{ii} = wavenumber_Guo(periods(ii),water_depth);

%We calculate L as shown in the manual and with the wavenumber calculated 
L{ii} = ((g*power(periods(ii),2))/(2*pi))*tanh(wavenumber{ii}.*water_depth); 
%ratio h/L
ratio_Lh{ii} = water_depth./L{ii};
%phase velocity with our function 
c{ii} = phase_velocity(L{ii},periods(ii)); 
%Propagation factor with our function 
prop_factor{ii} = propagation_factor(wavenumber{ii},water_depth); 
%group velocity with our function  
cg{ii} = group_velocity(c{ii},prop_factor{ii}); 
end


%Plot evolution of wavelength L as a function of depth 
figure() 
plot(water_depth, L{1}); 
title('Evolution of the wavelength L as a function of depth'); 
hold on 
plot(water_depth, L{2}); 
plot(water_depth, L{3}); 
grid on 
legend('T = 6 s', 'T = 9 s','T = 12 s'); 
xlabel('Water depth [m]') 
ylabel('Wavelength [m]') 
xlim([0 140])



% Plot evolution of the ratio h/L as a function of the depth

figure()
plot(water_depth,ratio_Lh{1}); 
xlim([0 140]); 
hold on
plot(water_depth,ratio_Lh{2}); 
plot(water_depth,ratio_Lh{3}); 
legend('T = 6 s', 'T = 9 s','T = 12 s'); 
title('Ratio h/L as a function of depth') 
grid on 
xlabel('Water depth [m]') 
ylabel('Ratio h/L') 
xlim([0 140])
%Adding horizontal lines with regions of water 
%Shallow water h/L <= 0.05 
h1 = yline(0.05);
%Deep water >= 0.5
h2 = yline(0.5);
set(get(get(h1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
text(80,0.03, 'Shallow water'); 
text(80,0.15, 'Intermediate water'); 
text(80,0.7, 'Deep water'); 

%Calculating the theoretical speed 
c_theory = sqrt(g.*water_depth);

%Evolution of C and Cg as a function of h 
figure()
subplot(2,1,1)
%Phase velocity 
plot(water_depth,c{1});
title('Change of phase velocity over depth with different periods')
hold on
plot(water_depth,c{2});
plot(water_depth,c{3});
%Theoretical c 
plot(water_depth,c_theory);
legend('T = 6 s', 'T = 9 s','T = 12 s','c = sqrt(gh)');
grid on 
xlabel('Water depth [m]')
ylabel('C [m/s]')
subplot(2,1,2)
%Group velocity 
plot(water_depth,cg{1});
title('Change of group velocity over depth with different periods')
hold on
plot(water_depth,cg{2});
plot(water_depth,cg{3});
%Theoretical c 
plot(water_depth,c_theory);
legend('T = 6 s', 'T = 9 s','T = 12 s','c = sqrt(gh)');
grid on 
xlabel('Water depth [m]')
ylabel('Cg [m/s]')

%% 3.2 Egmond Data 

%Loading data
mean_waterdepth = load('MeanWaterDepth.txt'); 
mean_waterdepth = mean_waterdepth';

%Periods from Table 9.1 of the manual 
%In order: low tide, mid tide and high tide 
periods_tides = [7.58 6.69 5.54]; %seconds 

%We run the loop for each tide 
for ii = 1:length(periods_tides(1,:))
    %We run the loop through each position 
      for jj = 1:length(mean_waterdepth(1,:))

%We calculate the wave number and L for each tide for all cross-shore
%positions with their correspondant period
wavenumber_Egmond(ii,jj) = wavenumber_Guo(periods_tides(ii),mean_waterdepth(ii,jj));
L_Egmond(ii,jj) = ((g*power(periods_tides(ii),2))/(2*pi))*tanh(wavenumber_Egmond(ii,jj).*mean_waterdepth(ii,jj));                       
    end
    
end

%Plot of wavelength for each cross-shore position for all tides
figure()
plot(L_Egmond(1,:),'-*'); 
hold on
plot(L_Egmond(2,:),'-*'); 
plot(L_Egmond(3,:),'-*'); 
xlabel('Cross-shore position')
xticklabels({'P1','P3','P4','P5','P6'}); 
ylabel('Wavelength [m]'); 
xticks([1 2 3 4 5])
legend('Low tide (T = 7.58 s)', 'Mid tide (T = 6.69 s)', 'High tide (T = 5.54 s)'); 
grid on

%% h/L ratio for Egmond data

%Calculation of ratios
ratio1=mean_waterdepth./L_Egmond(1,:);  
ratio2=mean_waterdepth./L_Egmond(2,:);
ratio3=mean_waterdepth./L_Egmond(3,:);

%Scatter plot with ratios for each period (and therefore tide) as a function of water depth. 
figure (1)
scatter(mean_waterdepth(1,:),ratio1(1,:));
hold on 
scatter(mean_waterdepth(2,:),ratio2(2,:))
hold on 
scatter(mean_waterdepth(3,:),ratio3(3,:))
%Adding horizontal lines with regions of water 
%Shallow water h/L <= 0.05 
yline(0.05, '--')
%Deep water >= 0.5
yline(0.5, '--')
xlabel('Water depth [m]')
ylabel('Ratio h/L')
ylim([0,0.55])
xlim([1, 5])
legend('T = 7.58 s (low tide)', 'T = 6.69 s (mid tide)', 'T = 5.54 s (high tide)','location','east')
title('Ratio h/L as a function of water depth for different periods T')
text(2.75,0.035,'Shallow')
text(2.75,0.065,'Intermediate')
text(2.75,0.49,'Intermediate')
text(2.75,0.515,'Deep')
grid on


