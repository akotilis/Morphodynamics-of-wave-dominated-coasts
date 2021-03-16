%% Chapter 7 

close all;
clear all;




%Compute and display the variance density spectra of the free-surface elevation for the three conditions. Choose the block size such that you have about 30 degrees of freedom. What do you observe?


nfft = length(condition1.data(:,1));

array_blocks = round(2*nfft/(12+1));

Fs = condition1.Fs; 

for ii = 1:3
    

[spectrum frequency edff confidenceint] = VarianceDensitySpectrum(eval(sprintf('condition%d.data(:,1)', ii)),array_blocks,Fs);

max_spectrum = max(spectrum);

spectrum_all{ii} = spectrum; 

confidence_intervals1{ii} = confidenceint(1)*spectrum; 
confidence_intervals2{ii} = confidenceint(2)*spectrum; 

%Sea swell wave height 
f_1 = 0.05; %Hz 
f_2 = Fs/2; %Nyquist frequency (sampling frequency/2) 
moment_swell = spectral_moment(frequency,spectrum,f_1,f_2,p)
Hmo_swell(ii) = 4*sqrt(moment_swell);

disp(sprintf('The sea swell wave height for condition %d is %d', ii, Hmo_swell(ii)))

%Infragravity wave height 
p= 0 ; %zero order m0 
f_1 = 0.005; %Hz 
f_2 = 0.05; %Hz 
moment_infra = spectral_moment(frequency,spectrum,f_1,f_2,p)
Hmo_inf(ii) = 4*sqrt(moment_infra);

disp(sprintf('The infragravity wave height for condition %d is %d', ii, Hmo_inf(ii)))

%Calculation of relative wave height 

H_relative(ii) = Hmo_swell(ii)./eval(sprintf('condition%d.data(:,1)', ii);

%Ratio infragravity over swell 

ratio = Hmo_inf(ii)/Hmo_swell(ii); 

%Removing the low frequency variations 
%filtered_data = fft_filter(eval(sprintf('condition%d.data(:,1)', ii), Fs, 0, fhigh)






figure(1)
subplot(3,1,ii)
plot(spectrum_all{1,ii})
title(sprintf('Variance density spectrum for condition %d ',ii))
hold on 
plot(confidence_intervals1{1,ii})
plot(confidence_intervals2{1,ii})
grid on 
xlabel('Frequency [Hz]'); 
ylabel('Variance spectral density [m^2/Hz]')
legend('Spectrum', 'Lower confidence interval', 'Upper confidence interval');

end


%% 7.2 


for ii = 1:3 
    
%Mean component of velocity u 
    
mean_u(ii) = mean(eval(sprintf('condition%d.data(:,2)', ii)));

%



    
    
    
end




%% 7.2 Sediment transport 


% 7.2.1 Preliminary analysis 

figure() 
subplot(1,2,1) 




%Sediment transport calculations 











