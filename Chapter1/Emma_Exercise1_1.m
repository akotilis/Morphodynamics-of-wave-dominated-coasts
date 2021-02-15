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

%We get the duration so we can make the time vector of our series 
duration_s = length/fs; %seconds
duration = (length/fs)/60; %min

%Calculating water depth h according to manual 

p = a1.*data+b1; %pressure [Pa]
h_a = p./(rho*g); %height above sensor [m]
h = h_a + h_b; %total depth water column [m]

h_mean = mean(h); %average water depth at sensor [m]

fprintf("The average water depth is %f\n meters", h_mean)

time = linspace(0.25,duration_s,length);



%Plotting h

figure()
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

%Detendring the data 
h_detrend = detrend(h);
%Removing the mean 
H = h - h_mean;

%Plotting the data without the mean vs detrend data 
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


