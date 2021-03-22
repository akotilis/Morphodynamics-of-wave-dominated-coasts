function [R beta] = velocity_acceleration_skewness(u,t)

%Function that calculates the degree of velocity skewness and acceleration
%skewness from the orbital velocity time series. 


%INPUT
%u --> orbital velocity (m/s) 
%t --> time (s) 

%OUTPUT 
%R --> velocity skewness (m/s) 
%beta --> acceleration skewness (m/sˆ2) 


u_crest = max(u); 

u_trough = min(u); 

R = u_crest/(u_crest + u_trough); 


end