%% Chapter 7 

close all;
clear all;


%% 7.1.1 Wave characteristics 

%Compute and display the variance density spectra of the free-surface elevation for the three conditions.
%Choose the block size such that you have about 30 degrees of freedom. What do you observe?

%We load the data file with the conditions and data 
load('ucData.mat')


nfft = length(condition1.data(:,1));

array_blocks = round(2*nfft/(12+1));

Fs = condition1.Fs; 

n = 3; %Number of conditions 

H_relative = zeros(nfft,3); 

for ii = 1:n
    
condition = eval(sprintf('condition%d.data(:,1)', ii)); 
    

[spectrum frequency edff confidenceint] = VarianceDensitySpectrum(condition,array_blocks,Fs);

max_spectrum = max(spectrum);

spectrum_all{ii} = spectrum; 

confidence_intervals1{ii} = confidenceint(1)*spectrum; 
confidence_intervals2{ii} = confidenceint(2)*spectrum; 

%Sea swell wave height 
p= 0 ; %zero order m0 
f_1 = 0.05; %Hz 
f_2 = Fs/2; %Nyquist frequency (sampling frequency/2) 
moment_swell = spectral_moment(frequency,spectrum,f_1,f_2,p)
Hmo_swell(ii) = 4*sqrt(moment_swell);

disp(sprintf('The sea swell wave height for condition %d is %d', ii, Hmo_swell(ii)))

%Infragravity wave height 
f_1 = 0.005; %Hz 
f_2 = 0.05; %Hz 
moment_infra = spectral_moment(frequency,spectrum,f_1,f_2,p)
Hmo_inf(ii) = 4*sqrt(moment_infra);

disp(sprintf('The infragravity wave height for condition %d is %d', ii, Hmo_inf(ii)))

%Calculation of relative wave height 

H_relative(ii) = Hmo_swell(ii)./eval(sprintf('condition%d.h', ii)); 

disp(sprintf('The relative wave height for condition %d is %d', ii, H_relative(ii)))


%Ratio infragravity over swell 

ratio(ii) = Hmo_inf(ii)/Hmo_swell(ii); 

disp(sprintf('The ratio between the infragavity wave height and sea swell wave height for condition %d is %d', ii, ratio(ii)))


%Removing the low frequency variations 
Fn = Fs/2; %Nyquist frequency
Flow = 0.05; %Hz 
filtered_data = fft_filter(condition, Fs, Flow, Fn); 


%Skewness and assymetry 

%[skewness assymetry] = skewness_asymmetry(filtered_data);


%Variance spectral density for all conditions 
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

%Skewness and assymetry for each condition 
figure(2)
subplot(3,1,ii) 
%plot(skewness)
hold on 
%plot(assymetry) 
title(sprintf('Skewness and assymetry for sea swell waves for condition %d',ii))
%legend('Skewness', 'Assymetry') 



end

%% 7.1.2 Separation of time series 

%Compute the mean component of the velocity, u , for the three conditions. 
%Is it consistent with the conclusions drawn in Section 7.1.1?

for ii = 1:n 
    
    u = eval(sprintf('condition%d.data(:,2)', ii)); 
    
    umean(ii) = mean(u); %cross-shore near bed velocity (m/s)
    
    disp(sprintf('The value of u mean is %0.2d for condition %d', umean(ii),ii)) 
    
    %Compute now the high and low frequency components for each condition using filter fft
    
    %u_hf = fft_filter(u,Fs, 0, ); 
    
    %u_lf = fft_filter(u,Fs,) 
    
    %Verifying that u = umean + u_hf + u_lf 
    
    %u_all = umean(ii) + u_hf(ii) + u_lf(ii); 
    
    %Comparison of u and umean + u_hf + u_lf 
    figure() 
    %plot(u)
    hold on 
    %plot(u_all) 
    title(sprintf('Comparison of near bed velocity (u) with umean + u_hf + u_lf for condition %d',ii))
    %legend('u', 'umean + u_hf + u_lf')
    
    
    %Figure of u_hf and u_lw for each condition 
    figure(5)
    subplot(3,1,ii) 
    %plot(u_hf) 
    hold on 
    %plot(u_lf) 
    %legend('Sea swell', 'Infragravity') 
    
end



%% 7.2 Sediment transport 


% 7.2.1 Preliminary analysis 

for ii = 1:n
    
%Velocity u - umean 

u = eval(sprintf('condition%d.data(:,2)', ii)); %cross-shore near bed velocity (m/s)

umean(ii) = mean(u); %Mean velocity (m/s) 

uminusumean = u - umean; % u - umean 

%u_lf = fft_filter(u,Fs,) 

concentration = eval(sprintf('condition%d.data(:,3)', ii)); %near-bed sediment concentration (kg/mˆ3)


figure() 
subplot(2,1,1) 
plot(uminusumean) 
%titulo = sprintf('Time evolution of u - $\\bar{u}$ for condition %d', ii)
%title(titulo,'Interpreter', 'LaTeX')
title(sprintf('Time evolution of u - u mean for condition %d', ii))
hold on
grid on 
ylabel('Cross-shore near bed velocity (m/s)') 
xlabel('Time (s)') 
xlim([0 4100])
%plot(u_lf) 
%legend('u - u mean', 'u low frequency') 
subplot(2,1,2) 
plot(concentration)
title(sprintf('Time evolution of the concentration for condition %d',ii)) 
grid on 
xlabel('Time (s)') 
ylabel('Concentration (kg/mˆ3)')
xlim([0 4100])

end


%Sediment transport calculations 

%Calculate the total sediment transport q as uc, i.e. the mean of 
%the product of the timeseries of u and c.

for ii = 1:n
    
 u = eval(sprintf('condition%d.data(:,2)', ii)) %cross-shore near bed velocity (m/s) 
 c = eval(sprintf('condition%d.data(:,3)', ii)) %near-bed sediment concentration (kg/mˆ3)

 %Total sediment transport 
total_sd(ii) = mean(u.*c); 

%Calculate the three parts of the sediment transport for each condition

%First term 
uc(ii) = mean(u)*mean(c); 


%low frequency velocity 
%u_hf = fft_filter(u,Fs, 0, ); 

%high frequency velocity  
%u_lf = fft_filter(u,Fs,) 

%low frequency concentration 
%c_hf = fft_filter(c,Fs, 0, ); 

%high frequency concentration  
%c_lf = fft_filter(c,Fs,) 

%Second term

%uhf_chf(ii) = mean(u_hf*c_hf); 

%Third term 
%ulf_clf = mean(u_lf*c_lf); 

%Adding the three terms 
all_terms(ii) = uc(ii) + uhf_chf(ii) + ulf_clf(ii); 

%Check that uc≈uc+ uhfchf + ulfclf.

figure()
subplot(3,1,ii)
plot(uc(ii))
hold on 
plot(all_terms(ii))
title(sprintf('Comparison ot total sediment transport and uc+ uhfchf + ulfclf for condition %d', ii)) 
legend('Total sediment transport','uc+ uhfchf + ulfclf')  


end











