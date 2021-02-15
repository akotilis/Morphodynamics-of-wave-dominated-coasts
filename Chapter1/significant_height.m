function [final_significant_wave_height] = significant_height(wave_heights)

%function [final_significant_wave_height] = significant_height(wave_heights)
% Calculation of significant wave height 
% input   wave_heights: input array of wave heights in m
% 
% output  final_significant_wave_height: array containing the value of significant
%         wave height of the wave_height series (sorting them on an
%         ascendent form and using H1/3 = 1/1*3n sum i = 2/3n+1 to n H_i

wave_heights = [1 1 2 3 3]; 
n = length(wave_heights); 
factor = 1/(n/3); 
ii = round((2*n/3) + 1); 
sorted_wave_heights = sort(wave_heights); 
significant_wave_height = 0; 
for jj = ii:n
significant_wave_height = significant_wave_height + sorted_wave_heights(jj);
end

final_significant_wave_height = factor*significant_wave_height; 
