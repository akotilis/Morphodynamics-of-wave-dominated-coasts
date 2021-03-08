%% Chapter 6 



% 6.1.2 

% Vector of Ur varying from 0.01 to 100

Ur = [0.01:0.01:100]; 

[skewness_R assymetry_R] = skewness_assymetryRuessink(Ur); 

figure()
subplot(2,1,1) 
semilogx(Ur,skewness_R) 
title('Skewness as a function of the Ursell number')
xlabel('Ursell number (log scale)')
ylabel('Skewness') 
grid on 

subplot(2,1,2) 
semilogx(Ur,assymetry_R) 
title('Assymetry as a function of the Ursell number')
xlabel('Ursell number (log scale)')
ylabel('Assymetry') 
grid on 

%% 6.2 Modelling Sk and As 





