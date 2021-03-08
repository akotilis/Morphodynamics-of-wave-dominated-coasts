function [sk as] = skewness_asymmetry(eta); 

%skewness_assymetry calculates the skewness and the assymetry from the 
%sea surface elevation as shown in Equations 6.1 and 6.2 of the manual.

%INPUT 
%eta ---> detrended sea surface elevation (meters) 

%OUTPUT 
%sk ---> skewness
%as ---> assymetry 

for ii = 1 :length(eta(1,:))
    
numerator = mean(eta(:,ii).^3); 
denominator = mean(eta(:,ii).^2).^(3/2); 

skewness(ii) = numerator/denominator; 

num_assymetry = mean(imag(hilbert(eta(:,ii)).^3)); 

assymetry(ii) = num_assymetry/denominator; 

sk(ii) = skewness(ii); 

as(ii) = assymetry(ii); 
end

end