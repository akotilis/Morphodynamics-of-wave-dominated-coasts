%% Chapter 8 

%% 8.1.1 Preliminary computations 

clear all; 
close all; 



T = 6; % period (seconds) 
Uw = 1 %orbital velocity amplitude (m/s)

r = [0 0.6]; 
phi = [0 -pi/2 -pi/4]; 

phi_names = {'0' '- \pi/2' '- \pi/4'}; 

%First case: r = 0 and phi = 0


[u{1} t{1}] = waveshape(r(1),phi(1),Uw,T); 
[R(1) beta(1) acceleration{1}] = velocity_acceleration_skewness(u{:,1},t{:,1}); 

figure(1) 
%subplot(4,1,1) 
plot(t{:,1},u{:,1}) 
title('Orbital velocity for different values of r  and \phi') 
%title(sprintf('Orbital velocity with r = %d and \\phi = %s ', round(r(1)),phi_names{1})); 
grid on 
%ylim([-1 2]) 
ylabel('u (m/s)') 
xlabel('Time (s)')

for ii = 1:length(phi)
    
%Calculating orbital velocity     

[u{ii+1} t{ii+1}] = waveshape(r(2),phi(ii),Uw,T); 

%Plot orbital velocity 
figure(1)
hold on; 
%subplot(4,1,ii+1) 
plot(t{:,ii+1},u{:,ii+1}) 
%title(sprintf('Orbital velocity with r = %.1d and \\phi = %s ',round(r(2)), phi_names{ii})); 
grid on 
ylim([-2 2])
%ylabel('u (m/s)') 
%xlabel('Time (s)')

%Calculation of R, beta and acceleration for all cases 
[R(ii+1) beta(ii+1) acceleration{ii+1}] = velocity_acceleration_skewness(u{:,ii+1},t{:,ii+1}); 



end

figure(1)
legend('r = 0 and \phi = 0', 'r = 0.6 and \phi = 0', 'r = 0.6 and \phi = - \pi/2','r = 0.6 and \phi = - \pi/4')


%Plots velocity and acceleration for cases 2 and 3 

%Case 2 
figure() 
subplot(2,1,1) 
plot(t{2}, u{2})
title('Time series of orbital velocity for r = 0.6 and \phi = 0 (case 2)') 
grid on 
ylim([-2 2]) 
ylabel('u (m/s)') 
xlabel('Time (s)')
legend('Velocity')
subplot(2,1,2) 
plot(t{2}, acceleration{2},'k') 
title('Time series of acceleration for r = 0.6 and \phi = 0 (case 2)') 
grid on 
ylim([-2.5 2.5]) 
ylabel('a (m/sˆ2)') 
xlabel('Time (s)')
legend('Acceleration')

%Case 3 
figure() 
subplot(2,1,1) 
plot(t{3}, u{3})
title('Time series of orbital velocity for r = 0.6 and \phi = -\pi/2 (case 3)') 
grid on 
ylim([-2 2]) 
ylabel('u (m/s)') 
xlabel('Time (s)')
legend('Velocity')
subplot(2,1,2) 
plot(t{3}, acceleration{3},'k') 
title('Time series of acceleration for r = 0.6 and \phi = -\pi/2 (case 3)') 
grid on
ylim([-2 2]) 
ylabel('a (m/sˆ2)') 
xlabel('Time (s)')
legend('Acceleration')







%% 8.1.2 Sediment transport analysis 








