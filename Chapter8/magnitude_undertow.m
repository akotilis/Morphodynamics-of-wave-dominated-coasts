function Ux = magnitude_undertow(E, Er, c, h, rho, Hrms) 

%Function magnitude_undertow calculates the magnitude of the undertow Ux as
%shown in Equation 8.5 of the manual (as Ux = M/rho*h_trough) 

%INPUT 
%E --> Wave energy 
%Er --> Roller energy 
%c --> phase velocity (m/s) 
%h --> Total water depth (meters) 
%rho --> density of water (kg/mˆ3) 
%Hrms --> Root mean square height (meters) 

%OUTPUT 
%Ux --> Magnitude of undertow 

M = E/c + 2*Er/c; 

rho = 1000 %kg/mˆ3 

h_trough = h - Hrms/2; 

Ux = M/(rho*h_trough); 



end 


