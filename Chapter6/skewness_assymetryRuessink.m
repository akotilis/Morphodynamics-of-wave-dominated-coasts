function [sk as] = skewness_assymetryRuessink(ursell_number)

% skewness_assymetryRuessink(ursell_number) computes the skewness and 
%assymetry from the Ruessink empirical fits defined in Equations 6.5 and
%6.6 of the manual. 

%Non linearity parameters 

inside_exp = (-0.471 - log10(ursell_number))./0.297; 

B = 0.857./(1 + exp(inside_exp)); 

Psi = -90 +90.*tanh(0.815./(ursell_number.^0.672)); 

%assymetry 

As = B.*sin(Psi); 

%skewness 

Sk = B.*cos(Psi);

sk = Sk; 

as = As; 



end