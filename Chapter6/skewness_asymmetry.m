function [sk as] = skewness_asymmetry(eta);



%skewness_assymetry calculates the skewness and the assymetry from the
%sea surface elevation as shown in Equations 6.1 and 6.2 of the manual.



%INPUT
%eta ---> detrended sea surface elevation (meters)



%OUTPUT
%sk ---> skewness
%as ---> assymetry




sk = (mean(eta.^3))./((mean(eta.^2)).^(3/2));



as = mean(imag(hilbert(eta)).^3)/(mean(eta.^2)).^(3/2);



end