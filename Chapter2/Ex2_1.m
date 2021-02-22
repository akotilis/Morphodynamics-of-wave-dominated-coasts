close all;
clear all;

%load data

low_tide = load('lowTide.txt');
mid_tide = load('midTide.txt');
high_tide = load('highTide.txt');


Fs = 2; %Sampling frequency Fs in Hz

 

%Calculating the density spectrum for last location on low tide using only
%one block

 %Defining one block

nfft = length(low_tide(:,1));

[S f edf conf95Interval] = VarianceDensitySpectrum(low_tide(:,1),nfft,Fs);


%Plotting the spectrum
figure(1) ; plot(f,S)
xlabel('frequency [Hz]')
ylabel('Variance spectral density [m^2/Hz]')
hold on
plot(f,conf95Interval(:,1)*S, 'red')
plot(f,conf95Interval(:,2)*S)
grid off ; legend('Variance density spectrum', 'lower 95% confidence interval', 'upper 95% confidence interval');
title('Variance spectral density for P1 at low tide and confidence intervals');

nfft = 2*length(low_tide(:,1))/(3+1);

%%Egmond data

nfft15 = length(low_tide(:,1))/8;

fE = zeros(5,2049);
SE = zeros(5,2049);
Ps = [1 3 4 5 6];
for ii = 1:5
    
    [S, f] = VarianceDensitySpectrum(low_tide(:,ii), nfft15, Fs);
    
    fE(ii, 1:length(f)) = f;
    SE(ii, 1:length(S)) = S;
   
    figure(2);  
    subplot(1,5,ii); plot(fE(ii,:), SE(ii,:));
    xlim([0 1]); ylim([0 2]);
    xlabel('Frequency(Hz)'); ylabel('Variance spectral density [m^2/Hz]');
    tit = sprintf('Low tide P %d', Ps(:,ii));
    title(tit);
    leg = sprintf('P%d', Ps(:,ii));
    legend(leg)
end


[S,f] = VarianceDensitySpectrum(low_tide(:,1), nfft15, Fs);

%Figure to estimate sea-swell frequency 
figure(3); plot(f,S);
xlim([0 1]); ylim([0 2]);