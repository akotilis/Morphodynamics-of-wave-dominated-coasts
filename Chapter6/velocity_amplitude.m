function Uw = velocity_amplitude(h,Hrms,T) 

%velocity_amplitude(h,Hrms,T) computes the velocity amplitude as Uw = pi/T
%* Hrms/sin(kh)

%INPUT 
%h ---> 
%Hrms ---> Root mean square wave height (meters) 
%T ---> Period (seconds) 

%OUTPUT 
%Uw ---> Velocity amplitude (meters) 

k = wavenumber_Guo(T,h); 

first = pi/T; 
second = Hrms/sinh(k*h); 

Uw = first*second; 

end