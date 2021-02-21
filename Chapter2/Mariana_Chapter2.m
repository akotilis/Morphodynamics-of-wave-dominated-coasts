%% Chapter 2 

close all; 
clear all; 

%Path where the data is located 
pathdata = "/Users/maribotton/Desktop/Wave_Coasts/morphodynamics_wavecoasts/Chapter1/Chapter1_Team7CassinoBeach/";

%Computation of the density spectrum with low tide 
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


% Plot on a new figure (vertical subplots) the spectra at the positions P1 to P6 
%for the low tide Egmond data. Ensure that all axes have the same range on the 
%x- and y-directions.

%% Part 2.2 Computation of spectral wave characteristics 

%Compute the total wave height Hm0 for the 5 time-series during all 3 tides
%(use the function spectral moment). 

%We will be using from now on 15 blocks 
nfft_15blocks = 2*nfft/(15+1);

array_names = ["lowTide.txt", "midTide.txt", "highTide.txt"]; 
n = length(array_names); %Number of datasets 
Fs = 2; %sampling frequency Hz 

%The first for loop runs through the datasets 
for jj = 1:n
    
    name = array_names(jj);
    wave_data = load(pathdata+name); 
    %Then we run the spectrum and moment for each data set for all
    %locations
for ii = 1:length(wave_data(1,:)) 
    
[spectrum frequency edff confidenceint] = VarianceDensitySpectrum(wave_data(:,ii),nfft_15blocks,Fs);

%Max spectrum to calculate peak period 
max_spectrum = max(spectrum); 

%We calculate the peak period for all spectrums and we save it 

peak_period = 1/max_spectrum; 

peak_period_all(ii,jj) = peak_period;  

%Total wave height 
p= 0 ; %zero order m0 
%For total wave height f1 = 0 and fn Nyquist (sampling frequency/2 i.e.
%Fs/2) 
f_1 = 0; 
f_2 = Fs/2; %Nyquist frequency (sampling frequency/2) 
moment = spectral_moment(frequency,spectrum,f_1,f_2,p);
Hm0(ii,jj) = 4*sqrt(moment); 

%Infragravity wave height 
f_1 = 0.005; %Hz 
f_2 = 0.05; %Hz 
moment_infra = spectral_moment(frequency,spectrum,f_1,f_2,p)
Hmo_inf(ii,jj) = 4*sqrt(moment_infra);

%Sea swell wave height 
f_1 = 0.05; %Hz 
f_2 = Fs/2; %Nyquist frequency (sampling frequency/2) 
moment_swell = spectral_moment(frequency,spectrum,f_1,f_2,p)
Hmo_swell(ii,jj) = 4*sqrt(moment_swell);

end

end

%Comparison between spectral wave height and significant wave height from
%last practical 


load(pathdata+'StatisticsEgmond.mat'); 



figure()
plot(Hm0(:,1), '-*','color', 'blue'); 
title('Spectral wave height and significant wave height for different cross-shore positions in Egmond') 
hold on 
plot(Hm0(:,2),'-*','color', 'red');
plot(Hm0(:,3),'-*','color', 'black');
%Significant wave height from chapter 1
plot(H13_total(:,1), '-o','color', 'blue'); 
plot(H13_total(:,2), '-o','color', 'red'); 
plot(H13_total(:,3), '-o','color', 'black'); 
legend('Low tide Hm0', 'Mid tide Hm0', 'High tide Hm0','Low tide H1/3', 'Mid tide H1/3', 'High tide H1/3'); 
xlabel('Cross-shore Position'); 
xticks([1 2 3 4 5])
xticklabels({'P1','P3','P4','P5','P6'}); 
ylabel('h[m]'); 
grid on 


%Cross-shore evolution of low tide data for infragravity and sea-swell for low tide

figure() 
plot(Hmo_inf(:,1),'-*'); 
title('Infragravity and sea-swell wave height for the low tide at Egmond')
hold on 
plot(Hmo_swell(:,1),'-*'); 
legend('Infragravity wave height', 'Sea-swell wave height');
xlabel('Cross-shore Position'); 
xticks([1 2 3 4 5])
xticklabels({'P1','P3','P4','P5','P6'}); 
ylabel('h[m]'); 
grid on 

save('Chapter2_waveheights','Hm0', 'Hmo_inf', 'Hmo_swell')

%Bonus question peak periods 
figure()
plot(peak_period_all(:,1),'-*'); 
title('Peak periods for cross-shore positions at Egmond'); 
hold on 
ylabel('Period [1/Hz]'); 
plot(peak_period_all(:,2),'-*');
plot(peak_period_all(:,3),'-*');
xlabel('Cross-shore Position'); 
legend('Low tide', 'Mid tide', 'High tide'); 
xticks([1 2 3 4 5])
xticklabels({'P1','P3','P4','P5','P6'}); 
grid on 




















