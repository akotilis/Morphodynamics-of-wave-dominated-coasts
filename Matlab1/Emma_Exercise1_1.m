clear all
close all

%Constants

fs = 4 %Hz
a1 = 19.97370 %calibration coef Pa/mV
b1 = -49.95959 %calibration coef Pa
rho = 1025 %kg/m3
g = 9.81 %m/s^2
h_b = 1.45 % height sea bed [m]

% load data
data = load('CalibP1.txt'); %mV

length = length(data);
duration_s = length/fs; %seconds
duration = (length/fs)/60; %min

%Calculating water depth h

p = a1.*data+b1; %pressure [Pa]
h_a = p./(rho*g); %height above sensor [m]
h = h_a + h_b; %total depth water column [m]

h_mean = mean(h); %average water depth at sensor [m]

fprintf("The average water depth is %f\n meters", h_mean)

time = linspace(0.25,duration_s,length);
%time = [1:1:duration_s]; %time vector [seconds]

%Plotting h

figure()
%plot(time,h(1:4:14400))
plot(time,h)
title('Water column depth as a function of time')
xlabel('Time [s]')
ylabel('Water depth [m]')
%xlim([0 3600])
hold on 
plot([time(1) time(end)],[h_mean h_mean], '--k')
grid on 
legend('h(t)', 'h_{mean}')
hold off

%Removing tide variations

h_detrend = detrend(h);
H = h - h_mean;
figure()
plot(time,H)
hold on
plot(time,h_detrend)
title('Water level as a function of time')
xlabel('Time [s]')
ylabel('Water depth [m]')
xlim([0 3600])
grid on 
legend('h(t) - h_{mean}', 'h_{detrend}')
hold off


