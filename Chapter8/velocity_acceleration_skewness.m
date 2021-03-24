function [R beta acceleration] = velocity_acceleration_skewness(u,t)

%Function that calculates the degree of velocity skewness and acceleration
%skewness from the orbital velocity time series. 


%INPUT
%u --> orbital velocity (m/s) 
%t --> time (s) 

%OUTPUT 
%R --> velocity skewness (m/s) 
%beta --> acceleration skewness (m/sˆ2) 
%acceleration --> acceleration from du/dt (m/sˆ2) 

%We calculate the amplitude between u_delta = 0 from the crest and from the
%trough for both velocity and acceleration 
u_crest = abs(max(u)); 

u_trough = abs(min(u)); 

%Calculation of parameter R 

R = u_crest./(u_crest + u_trough); 

%Gradient of velocity 

du = gradient(u);

%Time variable 

dt = t(1,2) - t(1,1); 

acceleration = du./dt; 

a_crest = abs(max(acceleration)); 

a_trough = abs(min(acceleration)); 

beta = a_crest./(a_crest + a_trough); 


end