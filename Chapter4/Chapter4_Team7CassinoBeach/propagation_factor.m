function propagation = propagation_factor(k,h)

%Function that computes the propagation factor from the wave number k and
%the water depth h

%Inputs: 
%k --> wave number 
%h --> water depth [meters]

%Output: 
%propagation --> propagation factor n 

mult = 2*k.*h; 

sine_part = mult./(sinh(mult)); 

n = 0.5.*(1+ sine_part); 

propagation = n; 


end