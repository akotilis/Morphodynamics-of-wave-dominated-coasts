%% Chapter 5 

%We first put the initial conditions for each tide to run the BJ Model 

%Angle of wave incidence 
theta = [-36 39 36]; %degrees 

%Period 1/3 
T0 = [7.58 6.69 5.54]; %seconds 

%Significant wave height 
H13 = [1.70 2.25 1.69]; %meters 

%Calculated root mean square wave height (as calculated H1/3 =
%sqrt(2)*Hrms)
Hrms0 = H13./sqrt(2); %meters 

%Water level (m NAP) 
Zeta = [-0.45 0.09 0.91] %meters 

%Bed profile for Egmond data 
bed_profile = load('prof1018.txt'); 

%first column: positions (x), second column: water level (z) 
profile = [bed_profile(:,1) bed_profile(:,2)]; 


% Model parameter 
hmin = 0.2;     % Minimal water depth for computation
                % (we stop the computation when h<hmin)


%Constants


ka = 0.022 % Apparent bed roughness [meters]

nu =  0.5  %large-scale mixing coefficient [mˆ2/s]
 for ii = 1:length(theta) 
     
waves = BJmodelEmma(Hrms0(ii),T0(ii),Zeta(ii),theta(ii),profile,hmin);

%Still missing definition of dzetady, wcross, wlong,


%longshoreCurrent(profile,dzetady,wcross,wlong,waves.c,theta(ii),waves.Dr,waves.ht,waves.st,ka,nu)


     
     
 end
 