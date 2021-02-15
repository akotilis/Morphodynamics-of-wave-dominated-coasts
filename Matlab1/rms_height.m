function [rms_waves] = rms_height(wave_heights) 
%function [rms_waves] = rms_height(wave_heights) 
% Calculation of root mean square height of the waves 
% input   wave_heights: input array of wave heights in m
% 
% output  rms_waves: array containing value of root mean square height of the waves
% defined as H_rms = sqrt(1/n * sum i = 1 to n (H_i)^2) 

n = length(wave_heights); 
factor = 1/n; 

sum_waveheights = 0; 
for jj = 1:n
    sum_waveheights = sum_waveheights + wave_heights(jj)*wave_heights(jj); 
end

rms_waves = sqrt(factor*sum_waveheights); 



end