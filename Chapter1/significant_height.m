function [final_significant_wave_height] = significant_height(wave_heights)

%function [final_significant_wave_height] = significant_height(wave_heights)
% Calculation of significant wave height 
% input   wave_heights: input array of wave heights in m
% 
% output  final_significant_wave_height: array containing the value of significant
%         wave height of the wave_height series (sorting them on an
%         ascendent form and using H1/3 = 1/1*3n sum i = 2/3n+1 to n H_i

n = length(wave_heights); 
factor = 1/(n/3); 
%From where the sum will go; we round it to the nearest integer to be able
%to perform the for 
ii = round((2*n/3) + 1); 
%We sort the wave heights 
sorted_wave_heights = sort(wave_heights); 
%We start the sum at zero 
significant_wave_height = 0; 
%We do the sum 
for jj = ii:n
significant_wave_height = significant_wave_height + sorted_wave_heights(jj);
end

%We multiply sum by factor 
final_significant_wave_height = factor*significant_wave_height; 
