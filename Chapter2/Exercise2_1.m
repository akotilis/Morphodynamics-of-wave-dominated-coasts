clear all
close all

% load data
data_lowtide = load('lowTide.txt');

%Constants
fs = 2; %sampling frequency [Hz]
dt = 1/fs; %delta t [s]

P1 = data_lowtide(:,1); %P1 is the sensor furthest away from the coast
nfft = length(P1); %Is the lenght of the block is as long as the signal, there can only be one block

%Check wether nBlocks is 1
%overlap = 0.5;
%n=size(P1,1);
%nOverlap = overlap*nfft;
%nBlocks = fix((n-nOverlap)/(nfft-nOverlap));

%Variance density function
[S f edf conf95Interval] = VarianceDensitySpectrum(P1, nfft, fs)

%Plotting the spectrum
plot(f,S)
xlabel('frequency [Hz]')
ylabel('Variance spectral density [m^2/Hz]')
hold on
plot(f,conf95Interval(:,1)*S)% 'k--', 'linewidth', 2)
hold on
plot(f,conf95Interval(:,2)*S)% 'k--', 'linewidth', 2)
legend('Spectral density', '95% confidence interval', '95% confidence interval')

%Compare with spectra obtained using 3,7,15 and 31 blocks
nfft3 = length(P1)/2
nfft7 = length(P1)/4
nfft15 = length(P1)/8
nfft31 = length(P1)/16

array_blocks = [nfft3 nfft7 nfft15 nfft31];
array_blocks = round(array_blocks);

for ii = 1:length(array_blocks)
    s = (-1)^array_blocks(ii);
    if s == -1
        array_blocks(ii) = array_blocks(ii) + 1;
    end 
    
    [S f edf conf95Interval] = VarianceDensitySpectrum(P1,array_blocks(ii),fs);
    
    max_spectrum = max(S);
    
    figure()
    plot(f,S)
    hold on
    ylim([0 max_spectrum])
end

sensors = ["P1", "P3", "P4", "P5", "P6"];

for i = 1:5
    
    [S f edf conf95Interval] = VarianceDensitySpectrum(data_lowtide(:,i),nfft15,fs);
    
    subplot(1,5,i)
    plot(f,S)
%     hold on
%     plot(f,conf95Interval(1)*S)
%     hold on
%     plot(f,conf95Interval(2)*S)
%     legend('Spectral density', 'lower confidence interval', 'upper confidence interval')
    title(sensors(i))
    xlabel('Frequency [Hz]')
    ylabel('Variance spectral density [m^2/Hz]')
    ylim([0,2])
end 