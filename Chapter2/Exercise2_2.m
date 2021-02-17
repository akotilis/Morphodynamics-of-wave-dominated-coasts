clear all
close all

% load data
data_lowtide = load('lowTide.txt');
data_hightide = load('highTide.txt'); 
data_midtide = load('midTide.txt'); 
%data_bedprofile = load('prof1018.txt');

%Constants
fs = 2; %sampling frequency [Hz]
dt = 1/fs; %delta t [s]
f_min = 0 %minimum frequency [Hz]
nfft15 = length(data_lowtide(:,1))/8

data = [data_lowtide, data_midtide, data_hightide];
names = ["Low Tide", "Mid Tide","High Tide"];
sensors = ["P1","P2","P3","P4","P5"];
moments = zeros(3,5);
H0 = zeros(3,5);

for i=1:3
    for j=1:5
        [S, f,] = VarianceDensitySpectrum(data(i),nfft15,fs);
        moments(i,j) = spectral_moment(f,S,0,max(f),0);
        H0(i,j) = 4*sqrt(moments(i,j));
    end
end
