function ur = ursell_number(k,h,Hrms)

%ursell_number calculates the Ursell number defined as Ur = 3ak/4(kh)ˆ3

%INPUT 
%k ---> wavenumber 
%Hrms ---> Root mean square wave height (meters) 
%h ---> water depth (meters) 

%OUTPUT 
%ur ---> ursell number 


amplitude = 0.5*sqrt(2).*Hrms; 

ursell = 3*amplitude.*k./(4*(k.*h).^3); 

ur = ursell; 


end