%% Chapter 2 

%Returns density spectrum S of signal x computed after dividing the signal
%into nfft samples using 50 percent of overlapping between samples (blocks)
%f is the frequency vector, degree of freedom edf 


%[S f edf conf95Interval] = VarianceDensitySpectrum(x,nfft,Fs)

pathdata = "/Users/maribotton/Desktop/Wave_Coasts/morphodynamics_wavecoasts/Chapter1/";


low_tide = load(pathdata+'lowTide.txt');
Fs = 2; %Sampling frequency Fs in Hz 

%Calculating the density spectrum for last location on low tide using only
%one block 

%Defining one block 

nfft = length(low_tide(:,1)); 

[spectrum_lowtide_lastlocation freq edf conf95Interval] = VarianceDensitySpectrum(low_tide(:,1),nfft,Fs);

%Visualization of the spectrum with one block 

conf95Interval2(:,1) = conf95Interval(:,1)*spectrum_lowtide_lastlocation;
conf95Interval2(:,2) = conf95Interval(:,2)*spectrum_lowtide_lastlocation;

figure()
plot(spectrum_lowtide_lastlocation) 
title('Spectrum of frequencies for P1 for low tide in Egmond');
hold on 
plot(conf95Interval2(:,1)); 
plot(conf95Interval2(:,2)); 
grid on 
xlim([0 2050]) 
xlabel('Frequency [Hz]'); 
legend('Spectrum', 'Lower confidence interval', 'Upper confidence interval'); 
ylabel('Variance spectral density [m^2/Hz]')


%Now doing it with 3,7,15,31 blocks 

%3 Blocks 

array_blocks = [2*nfft/(3+1) 2*nfft/(7+1) 2*nfft/(15+1) 2*nfft/(31+1)];



for ii = 1:length(array_blocks) 

[spectrum frequency edff confidenceint] = VarianceDensitySpectrum(low_tide(:,1),array_blocks(ii),Fs);

max_spectrum = max(spectrum);

spectrum_all{ii} = spectrum; 

confidence_intervals1{ii} = confidenceint(1)*spectrum; 
confidence_intervals2{ii} = confidenceint(2)*spectrum; 

end

figure() 
%3 Blocks 
subplot(2,2,1) 
plot(spectrum_all{1,1})
title('3 Blocks')
hold on 
plot(confidence_intervals1{1,1})
plot(confidence_intervals2{1,1})
grid on 
xlabel('Frequency [Hz]'); 
ylabel('Variance spectral density [m^2/Hz]')
legend('Spectrum', 'Lower confidence interval', 'Lower confidence interval'); 
ylim([0 10]) 
%7 Blocks 
subplot(2,2,2) 
plot(spectrum_all{1,2})
title('7 Blocks')
hold on 
plot(confidence_intervals1{1,2})
plot(confidence_intervals2{1,2})
grid on 
xlabel('Frequency [Hz]'); 
ylabel('Variance spectral density [m^2/Hz]')
ylim([0 10]) 
legend('Spectrum', 'Lower confidence interval', 'Upper confidence interval'); 
%15 Blocks 
subplot(2,2,3) 
plot(spectrum_all{1,3})
title('15 Blocks')
hold on 
plot(confidence_intervals1{1,3})
plot(confidence_intervals2{1,3})
grid on 
xlabel('Frequency [Hz]'); 
ylabel('Variance spectral density [m^2/Hz]')
ylim([0 10]) 
legend('Spectrum', 'Lower confidence interval', 'Upper confidence interval'); 
%31 Blocks 
subplot(2,2,4)
plot(spectrum_all{1,4})
title('31 Blocks')
hold on 
plot(confidence_intervals1{1,4})
plot(confidence_intervals2{1,4})
grid on 
xlabel('Frequency [Hz]'); 
ylabel('Variance spectral density [m^2/Hz]')
ylim([0 10]) 
legend('Spectrum', 'Lower confidence interval', 'Upper confidence interval'); 


%% Part 2.2 Computation of spectral wave characteristics 

%Compute the total wave height Hm0 for the 5 time-series during all 3 tides (use the function spectral moment). 


%Store the results in a 5 × 3 array.
















