clear all;
close all;

% Definition of cross-shore coordinates (m)
x = (1:1:500)';  

% Definition of zb, bed elevation relative to mean water level (m)
zb = -9 * ones(size(x));   
zb = zb + (x<=260).*(x-160)/20  + ((x<=300)&(x>260)).*(-1/40*(x-260) + 5) + (x>300).*(1/20*(x-300) + 4); 

% Definition of the array profile, input argument for BJmodel
profile = [x zb];

% Offshore wave conditions
Hrms0 = 1;      % Root mean square wave height (m)
theta0 = 0;     % Angle of incidence (degrees)
T0 = 10;        % Characteristic period (s)
Zeta = 0;       % Mean water level (m)

% Model parameter 
hmin = 0.2;     % Minimal water depth for computation
                % (we stop the computation when h<hmin)

%----------------------------------------------------------------------------------
%       Computation and visualisation of wave characteristics using a variable zeta
%----------------------------------------------------------------------------------
%waves = BJmodelEmma(Hrms0,T0,Zeta,theta0,profile,hmin);

Zeta = [-1 0 1]
Title = ["Zeta = -1m", "Zeta = 0m", "Zeta = 1m"];

for i=1:length(Zeta);
    waves(i) = BJmodelEmma(Hrms0,T0,Zeta(i),theta0,profile,hmin);
    figure(i);
    subplot(5,1,1); 
    plot(waves(i).x,waves(i).Hrms)
    title(Title(i))
    ylabel('Hrms (m)')
    xlim([0,500]);
    subplot(5,1,2); plot(waves(i).x,waves(i).z,'k')
    hold on; plot(waves(i).x,Zeta(i)*ones(size(x)),'-.')
    ylabel('zb (m)')
    xlim([0,500]);
    subplot(5,1,3); plot(waves(i).x,waves(i).Dbr)
    ylabel('D_{Br} (W/m^2)')
    xlim([0,500]);
    subplot(5,1,4); plot(waves(i).x,waves(i).Dr)
    xlim([0,500]);
    ylabel('D_r (W/m^2)')
    subplot(5,1,5); plot(waves(i).x,waves(i).eta)
    xlim([0,500]);
    xlabel('x (m)')
    ylabel('eta (m)')
end 


%----------------------------------------------------------------------------------
%       Computation and visualisation of wave characteristics using a
%       variable Hrms
%----------------------------------------------------------------------------------


Zeta = 0
Hrms0 = [0.5, 2]
Title = ["Hrms0 = 0.5m", "Hrms0 = 2m"];

for i=1:length(Hrms0);
    waves(i) = BJmodelEmma(Hrms0(i),T0,Zeta,theta0,profile,hmin);
    figure(i);
    subplot(5,1,1); plot(waves(i).x,waves(i).Hrms)
    title(Title(i))
    ylabel('Hrms (m)')
    xlim([0,500]);
    subplot(5,1,2); plot(waves(i).x,waves(i).z,'k')
    hold on; plot(waves(i).x,Zeta*ones(size(x)),'-.')
    ylabel('zb (m)')
    xlim([0,500]);
    subplot(5,1,3); plot(waves(i).x,waves(i).Dbr)
    ylabel('D_{Br} (W/m^2)')
    xlim([0,500]);
    subplot(5,1,4); plot(waves(i).x,waves(i).Dr)
    xlim([0,500]);
    ylabel('D_r (W/m^2)')
    subplot(5,1,5); plot(waves(i).x,waves(i).eta)
    xlim([0,500]);
    xlabel('x (m)')
    ylabel('eta (m)')
end
