function propagation = propagation_factor(k,h)

%Function that computes the propagation factor from the wave number k and
%the water depth h

%Inputs: 
%k --> wave number 
%h --> water depth [meters]

%Output: 
%propagation --> propagation factor n 

sine_part = 2*k*h/(sinh(2*k*h)); 

n = 0.5*(1+ sine_part); 

propagation = n; 


end