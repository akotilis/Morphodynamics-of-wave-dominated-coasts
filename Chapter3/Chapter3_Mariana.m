%% Chapter 3 : Linear wave theory 

clear all; 
close all; 

periods = [6 9 12]; %seconds 
water_depth = [0:0.5:140]; %meters

%Constants 

g  = 9.81; %Gravity, m/sˆ2 
omega = 2*pi./periods;

for ii = 1:length(periods)
wavenumber{ii} = wavenumber_Guo(periods(ii),water_depth);
L{ii} = ((g*power(periods(ii),2))/(2*pi))*tanh(wavenumber{ii}.*water_depth); 
%ratio h/L
ratio_Lh{ii} = water_depth./L{ii};
%phase velocity
c{ii} = phase_velocity(L{ii},periods(ii)); 
%group velocity 
cg{ii} = group_velocity(omega(ii),wavenumber{ii}); 
end

%Plot evolution of wavelength L as a function of depth 
figure() 
plot(water_depth, L{1}); 
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
yline(0.05)
%Deep water >= 0.5
yline(0.5)


%Evolution of C and Cg as a function of h 
figure()
subplot(2,1,1)
%Phase velocity 
plot(water_depth,c{1});
title('Change of phase velocity over depth with different periods')
hold on
plot(water_depth,c{2});
plot(water_depth,c{3});
legend('T = 6 s', 'T = 9 s','T = 12 s');
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
legend('T = 6 s', 'T = 9 s','T = 12 s');
grid on 
xlabel('Water depth [m]')
ylabel('Cg [m/s]')


%% 3.2 Egmond Data 

mean_waterdepth = load('MeanWaterDepth.txt'); 
%Periods from Table 9.1 of the manual 

%In order: low tide, mid tide and high tide 

periods_tides = [7.58 6.69 5.54]; %seconds 

%We run the loop for each tide 
for ii = 1:length(mean_waterdepth(1,:))
    %We run the loop through each position 
      for jj = 1:length(mean_waterdepth(:,1))
          
wavenumber_Egmond(jj,ii) = wavenumber_Guo(periods_tides(ii),mean_waterdepth(jj,ii));
L_Egmond(jj,ii) = ((g*power(periods_tides(ii),2))/(2*pi))*tanh(wavenumber_Egmond(jj,ii).*mean_waterdepth(jj,ii)); 
                        
    end
    
end

figure()
plot(L_Egmond(:,1),'-*'); 
hold on
plot(L_Egmond(:,2),'-*'); 
plot(L_Egmond(:,3),'-*'); 
xlabel('Cross-shore position')
xticklabels({'P1','P3','P4','P5','P6'}); 
ylabel('Wavelength [m]'); 
xticks([1 2 3 4 5])
legend('Low tide (T = 7.58 s)', 'Mid tide (T = 6.69 s)', 'High tide (T = 5.54 s)'); 
grid on 

















