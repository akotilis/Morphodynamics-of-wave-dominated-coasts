function wavenumber = wavenumber_Guo(T,h) 

%Computes wave number with the approximation of Guo (2002). 
%Inputs: 
%T --> Wave period in seconds 
%h --> Water depth in meters 

%Output: 
%Wave number 

omega = 2*pi/T;
g = 9.81; %Gravity constant [m/sˆ2] 

x = h*omega./sqrt(g*h); 

beta = 2.4908; 

exponential = exp(-power(x,beta)); 

y = (x.*x).*power(1 - exponential,-1/beta);  

k = y./h;

wavenumber = k; 

end