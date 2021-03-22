%% Chapter 8 

%% 8.1.1 Preliminary computations 


T = 6; % period (seconds) 
Uw = 1 %orbital velocity amplitude (m/s)

r = [0 0.6]; 
phi = [0 -pi/2 -pi/4]; 

%First case: r = 0 and phi = 0

[u{1} t{1}] = waveshape(r(1),phi(1),Uw,T); 


for ii = 1:length(phi)

[u{ii+1} t{ii+1}] = waveshape(r(2),phi(ii),Uw,T); 

end

figure()
plot(t{1},u{1})

disp(max(u{1}))

%% 


