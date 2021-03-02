function wavenumber = wavenumber_Guo(T,h) 

%Computes wave number with the approximation of Guo (2002). 
%Inputs: 
%T --> Wave period in seconds 
%h --> Water depth in meters 

%Output: 
%Wave number 


omega = 2*pi/T; %Frequency (1/second) 
g = 9.81; %Gravity constant [m/sˆ2] 

x = h*omega./sqrt(g*h); %As shown in the manual 

beta = 2.4908; %value from the manual 

exponential = exp(-power(x,beta)); %We calculate the exponential term separately to facilitate 
%the calculation of y 

y = (x.*x).*power(1 - exponential,-1/beta);  

%We calculate the wavenumber
k = y./h;

%We return the wave number 
wavenumber = k; 

end