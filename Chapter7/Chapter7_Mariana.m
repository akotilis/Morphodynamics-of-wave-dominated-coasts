%% Chapter 7 

close all;
clear all;


%% 7.1.1 Wave characteristics 

%Compute and display the variance density spectra of the free-surface elevation for the three conditions.
%Choose the block size such that you have about 30 degrees of freedom. What do you observe?

%We load the data file with the conditions and data 
load('ucData.mat')


%We do a first overview of the detrended sea surface elevation for the 3
%sets of data 
for ii = 1:3 
figure(50)
subplot(3,1,ii)
%Plot free surface elevation 
plot(eval(sprintf('condition%d.data(:,1)', ii)))
grid on 
xlabel('Time (s)')
ylabel('z (m)')
ylim([-2.2 2.2]) 
xlim([0 4096])
title(sprintf('Detrended sea surface elevation for condition %d',ii)) 

end


%Defining length of the data 
nfft = length(condition1.data(:,1));

%We define 12 blocks to have 30 degrees of freedom 
array_blocks = round(2*nfft/(12+1));

%We define Fs (is 2 for all data sets) 
Fs = condition1.Fs; %Hz

n = 3; %Number of conditions 

%We run for all conditions (all data sets) for the first dimension (i.e.
%the detrended sea surface elevation 
for ii = 1:n
 
%We read the sea surface elevation data     
condition = eval(sprintf('condition%d.data(:,1)', ii)); 
    
%We calculate the variance spectrum 
[spectrum frequency edff confidenceint] = VarianceDensitySpectrum(condition,array_blocks,Fs);

%We save the spectrum 
spectrum_all{ii} = spectrum; 

%Sea swell wave height 
p= 0 ; %zero order m0 
f_1 = 0.05; %Hz 
f_2 = Fs/2; %Nyquist frequency (sampling frequency/2) 
%We calculate the spectral moment 
moment_swell = spectral_moment(frequency,spectrum,f_1,f_2,p)
%Calculation of sea swell wave height 
Hmo_swell(ii) = 4*sqrt(moment_swell);

disp(sprintf('The sea swell wave height for condition %d is %d', ii, Hmo_swell(ii)))

%Infragravity wave height 
f_1 = 0.005; %Hz 
f_2 = 0.05; %Hz 
%We calculate the spectral moment 
moment_infra = spectral_moment(frequency,spectrum,f_1,f_2,p)
%Calculation of infragravity wave height 
Hmo_inf(ii) = 4*sqrt(moment_infra);

disp(sprintf('The infragravity wave height for condition %d is %d', ii, Hmo_inf(ii)))

%Calculation of relative wave height 

H_relative(ii) = Hmo_swell(ii)/eval(sprintf('condition%d.h', ii)); 

disp(sprintf('The relative wave height for condition %d is %d', ii, H_relative(ii)))


%Ratio infragravity over swell 

ratio(ii) = Hmo_inf(ii)/Hmo_swell(ii); 

disp(sprintf('The ratio between the infragavity wave height and sea swell wave height for condition %d is %d', ii, ratio(ii)))


%Removing the low frequency variations 
Fn = Fs/2; %Nyquist frequency
Flow = 0.05; %Hz 
%We filter the data for the low frequencies 
filtered_data = fft_filter(condition, Fs, Flow, Fn); 


%Skewness and assymetry calculation with the filtered data  

[skewness(ii) assymetry(ii)] = skewness_asymmetry(filtered_data);


%Variance spectral density for all conditions 
figure(1)
subplot(3,1,ii)
plot(frequency,spectrum_all{1,ii})
title(sprintf('Variance density spectrum for condition %d ',ii))
grid on 
xlabel('Frequency [Hz]'); 
ylabel('Var. spec. density [m^2/Hz]')
xlim([0 1]) 
end

%% 7.1.2 Separation of time series 

%Compute the mean component of the velocity, u , for the three conditions. 
%Is it consistent with the conclusions drawn in Section 7.1.1?

%Running for all data sets (all conditions) 
for ii = 1:n 
    
    %We read the velocity from the data set 
    u = eval(sprintf('condition%d.data(:,2)', ii)); 
    
    %We calculate the mean of the velocity 
    umean(ii) = mean(u); %cross-shore near bed velocity (m/s)
    
    disp(sprintf('The value of u mean is %0.2d for condition %d', umean(ii),ii)) 
    
    %Compute now the high and low frequency components for each condition 
    %using filter fft
   
    u_hf(:,ii) = fft_filter(u,Fs, Flow,Fn); 
    
    u_lf(:,ii) = fft_filter(u,Fs,0,Flow); 
    
    %Verifying that u = umean + u_hf + u_lf 
    u_components(:,ii) = u_hf(:,ii) + u_lf(:,ii); 
    
    u_all(:,ii)= u_components(:,ii) + umean(ii);
    
    %Comparison of u and umean + u_hf + u_lf 
    figure(3)
    subplot(3,1,ii)
    plot(u)
    hold on 
    plot(u_all(:,ii)) 
    plot(u-u_all(:,ii))
    xlim([0 4096])
    title(sprintf('Comparison of near bed velocity (u) with umean + u_hf + u_lf for condition %d',ii))
    legend('u', 'u_m_e_a_n + u_h_f + u_l_f', 'difference')
    xlabel('Time (s)')
    ylabel('Velocity (m/s)') 
    grid on
    
    
    %Figure of u_hf and u_lf for each condition 
    figure(5)
    subplot(3,1,ii) 
    plot(u_hf(:,ii)) 
    hold on 
    plot(u_lf(:,ii)) 
    title(sprintf('Sea swell and infragravity velocity components for condition %d',ii))
    legend('Sea swell (u_h_f)', 'Infragravity (u_l_f)') 
    xlabel('Time (s)')
    ylabel('Velocity (m/s)') 
    xlim([0 4096]) 
    grid on
    
end



%% 7.2 Sediment transport 


% 7.2.1 Preliminary analysis 

%We run for all datasets 
for ii = 1:n
    
%Velocity u - umean 

%We read the velocity from the data 
u = eval(sprintf('condition%d.data(:,2)', ii)); %cross-shore near bed velocity (m/s)

%We calculate the mean of the velocity 
umean(ii) = mean(u); %Mean velocity (m/s) 

%We do u minus umean 
uminusumean = u - umean; % u - umean 

%We read the concentration from the data 
concentration = eval(sprintf('condition%d.data(:,3)', ii)); %near-bed sediment concentration (kg/m??3)

%Plot of the u minus umean and u_lf, u_hf (top) and concentration (bottom) 
figure() 
subplot(2,1,1) 
plot(uminusumean(:,ii)) 
title(sprintf('Time evolution of u - u_m_e_a_n for condition %d', ii))
hold on
grid on 
ylabel('Cross-shore near bed velocity (m/s)') 
xlabel('Time (s)') 
xlim([0 300])
plot(u_lf(:,ii),'k') 
plot(u_hf(:,ii))
legend('u - u_m_e_a_n', 'Infragravity velocity (u_l_f)', 'Sea swell velocity (u_h_f)') 
subplot(2,1,2) 
plot(concentration)
title(sprintf('Time evolution of the concentration for condition %d',ii)) 
grid on 
xlabel('Time (s)') 
ylabel('Concentration (kg/m??3)')
xlim([0 300])

end


%Sediment transport calculations 

%Calculate the total sediment transport q as uc, i.e. the mean of 
%the product of the timeseries of u and c.

%We run for all datasets 
for ii = 1:n
 
 %We read the velocity 
 u = eval(sprintf('condition%d.data(:,2)', ii)) %cross-shore near bed velocity (m/s)
 %We read the concentration 
 c = eval(sprintf('condition%d.data(:,3)', ii)) %near-bed sediment concentration (kg/m??3)

%Total sediment transport 
q(ii) = mean(u.*c); 

%Calculate the three parts of the sediment transport for each condition

%First term 
uc(ii) = mean(u).*mean(c); 

%high frequency concentration  
c_hf(:,ii) = fft_filter(c,Fs,Flow,Fn); 

%low frequency concentration 
c_lf(:,ii) = fft_filter(c,Fs,0,Flow); 

%Second term

uhf_chf(ii) = mean(u_hf(:,ii).*c_hf(:,ii)); 

%Third term 
ulf_clf(:,ii) = mean(u_lf(:,ii).*c_lf(:,ii)); 

%Adding the three terms 
all_terms(ii) = uc(ii) + uhf_chf(ii) + ulf_clf(ii); 

%Check that uc???uc+ uhfchf + ulfclf.

end











