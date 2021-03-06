clear all
close all

%Constants
fs = 2; %sampling frequency [Hz]

%Load data
low_tide = load('lowTide.txt'); %sea surface elevation [m]
high_tide = load('highTide.txt'); %sea surface elevation [m]

duration = length(low_tide)/fs; %seconds
length_data = length(high_tide(:,1)); 
time = linspace(0.25,duration,length_data);


figure()
subplot(3,2,1)
%Plotting the low tide for P1 position 
plot(time, low_tide(:,1))
xlabel('Time [s]')
ylabel('h [m]')
title('Low Tide P1')
grid on
xlim([0 time(end)])
ylim([-2 2])
%Plotting the low tide for P3 position 
subplot(3,2,3)
plot(time, low_tide(:,2))
xlabel('Time [s]')
ylabel('h [m]')
title('Low Tide P3')
grid on
xlim([0 time(end)])
ylim([-2 2])
%Plotting the low tide for P6 position 
subplot(3,2,5)
plot(time, low_tide(:,5))
xlabel('Time [s]')
ylabel('h [m]')
title('Low Tide P6')
grid on
xlim([0 time(end)])
ylim([-2 2])
%Plotting the high tide for P1 position 
subplot(3,2,2)
plot(time, high_tide(:,1))
xlabel('Time [s]')
ylabel('h [m]')
title('High Tide P1')
grid on
xlim([0 time(end)])
ylim([-2 2])
%Plotting the high tide for P3 position 
subplot(3,2,4)
plot(time, high_tide(:,2))
xlabel('Time [s]')
ylabel('h [m]')
title('High Tide P3')
grid on
xlim([0 time(end)])
ylim([-2 2])
%Plotting the high tide for P6 position 
subplot(3,2,6)
plot(time, high_tide(:,5))
xlabel('Time [s]')
ylabel('h [m]')
title('High Tide P6')
grid on
xlim([0 time(end)])
ylim([-2 2])





