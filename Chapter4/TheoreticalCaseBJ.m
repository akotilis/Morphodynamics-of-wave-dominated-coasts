% Main script for the computation of the cross-shore wave evolution 
% for a simplified bar-trough case based on the Battjes and Janssen (1978) model 

%------------------------------------
%           Initialisation
%------------------------------------

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

Zeta = [-1 0 1]; %meters
theta0 = 0;
Hrms0 = 1;

Title = ["Zeta = -1m", "Zeta = 0m", "Zeta = 1m"];

for i=1:length(Zeta);
    waves(i) = BJmodelEmma(Hrms0,T0,Zeta(i),theta0,profile,hmin);    
end 

% Figure for 4.2 (3rd task)
figure(); 
subplot(5,1,1); plot(waves(2).x, waves(2).Hrms); ylabel('Hrms(m)');
xlim([0 500]); ylim([0 2]); hold on;title('Water level \zeta = 0m');

subplot(5,1,2); plot(waves(2).x, waves(2).eta); ylabel('set-up (m)');xlim([0 500]); ylim([-0.2 0.2]); hold on;
subplot(5,1,3); plot(waves(2).x, waves(2).Dbr); ylabel('D_{br} (W/m^2)');xlim([0 500]); ylim([0 400]); hold on;
subplot(5,1,4); plot(waves(2).x, waves(2).Dr); ylabel('D_{r} (W/m^2)');xlim([0 500]); ylim([0 400]); hold on;
subplot(5,1,5); plot(waves(2).x, waves(2).z, 'k'); hold on; plot(waves(2).x, Zeta(2)*ones(size(x)), '-.'); ylabel('zb (m)'); 
xlim([0 500]); ylim([-20 10]); xlabel('x(m)'); 
%%
% Figure for 4.2(4th task)
figure();
subplot(5,1,1); plot(waves(1).x,waves(1).Hrms) 
hold on; plot(waves(2).x, waves(2).Hrms)
hold on; plot(waves(3).x, waves(3).Hrms) 
legend('\zeta = -1', '\zeta =0', '\zeta = 1'); legend('location','northeast');
title('Water level \zeta = -1m, \zeta=0m and \zeta= 1m')
ylabel('Hrms (m)')
xlim([0,500]); ylim([0 2]); 

subplot(5,1,2); plot(waves(1).x,waves(1).eta)
hold on; plot(waves(2).x, waves(2).eta)
hold on; plot(waves(3).x, waves(3).eta)
xlim([0 500]); ylim([-0.2 0.2]);
xlabel('x (m)'); ylabel('set-up (m)');



subplot(5,1,3); plot(waves(1).x,waves(1).Dbr)
hold on; plot(waves(2).x, waves(2).Dbr)
hold on; plot(waves(3).x, waves(3).Dbr)
ylabel('D_{Br} (W/m^2)')
xlim([0,500]); ylim([0 400]);

subplot(5,1,4); plot(waves(1).x,waves(1).Dr)
hold on; plot(waves(2).x, waves(2).Dr)
hold on; plot(waves(3).x, waves(3).Dr)
xlim([0,500]); ylim([0 400]);
ylabel('D_r (W/m^2)')

% subplot(5,1,5); plot(waves(1).x,waves(1).eta)
% hold on; plot(waves(2).x, waves(2).eta)
% hold on; plot(waves(3).x, waves(3).eta)

subplot(5,1,5); plot(waves(1).x,waves(1).z,'k')
hold on; plot(waves(2).x, waves(2).z, 'k')
hold on; plot(waves(3).x, waves(3).z, 'k')
hold on; plot(waves(1).x, Zeta(1)*ones(size(x)),'-.')
hold on; plot(waves(2).x, Zeta(2)*ones(size(x)), '-.')
hold on; plot(waves(3).x, Zeta(3)*ones(size(x)),'-.')
ylabel('zb (m)')
xlim([0,500]); ylim([-20 10]);

%%
%----------------------------------------------------------------------------------
%       Computation and visualisation of wave characteristics using a
%       variable Hrms
%----------------------------------------------------------------------------------


Zeta = 0 %meters
Hrms0 = [0.5, 2]; %meters
theta0=0;
Title = ["Hrms0 = 0.5m", "Hrms0 = 2m"];

for i=1:length(Hrms0)
    waves(i) = BJmodelEmma(Hrms0(i),T0,Zeta,theta0,profile,hmin);
end 

figure();
subplot(5,1,1); plot(waves(1).x,waves(1).Hrms)
hold on; plot(waves(2).x, waves(2).Hrms)
title('Hrms = 0.5m and 2 m, \zeta =0m')
ylabel('Hrms (m)'); legend('Hrms=0.5', 'Hrms=2')
xlim([0,500]); ylim([0 3])
subplot(5,1,5); plot(waves(1).x,waves(1).z,'k')
hold on; plot(waves(2).x, waves(2).z, 'k')
hold on; plot(waves(1).x,Zeta*ones(size(x)),'-.')
hold on; plot(waves(1).x,Zeta*ones(size(x)),'-.')
ylabel('zb (m)')
xlim([0,500]);
subplot(5,1,3); plot(waves(1).x,waves(1).Dbr)
hold on; plot(waves(2).x, waves(2).Dbr)
ylabel('D_{Br} (W/m^2)')
xlim([0,500]);
subplot(5,1,4); plot(waves(1).x,waves(1).Dr)
hold on; plot(waves(2).x, waves(2).Dr)
xlim([0,500]);
ylabel('D_r (W/m^2)')
subplot(5,1,2); plot(waves(1).x,waves(1).eta)
hold on; plot(waves(2).x, waves(2).eta)
xlim([0,500]); ylim([-0.2 0.4])
xlabel('x (m)')
ylabel('set-up (m)')


%----------------------------------------------------------------------------------
%       Computation and visualisation of wave characteristics using a
%       variable theta
%----------------------------------------------------------------------------------
%%
Zeta = 0 %meters
Hrms0 = 1 %meters
theta0 = [22.5, 45]; %degrees
Title = ["theta0 = 22.5 degrees", "theta0 = 45 degrees"];

for i=1:length(theta0);
    waves(i) = BJmodelEmma(Hrms0,T0,Zeta,theta0(i),profile,hmin);
end

figure();
subplot(5,1,1); plot(waves(1).x,waves(1).Hrms)
hold on; plot(waves(2).x, waves(2).Hrms)
title('\theta = 22.5^{\circ} and \theta=45^{\circ}') ; legend('\theta = 22.5^{\circ}', '\theta = 45^{\circ}');
ylabel('Hrms (m)')
xlim([0,500]); ylim([0 2]);

subplot(5,1,2); plot(waves(1).x,waves(1).eta)
hold on; plot(waves(2).x, waves(2).eta)
xlim([0,500]); ylim([-0.2 0.2])
ylabel('set-up (m)')

subplot(5,1,3); plot(waves(1).x,waves(1).Dbr)
hold on; plot(waves(2).x, waves(2).Dbr)
ylabel('D_{Br} (W/m^2)')
xlim([0,500]);ylim([0 400]);

subplot(5,1,4); plot(waves(1).x,waves(1).Dr)
hold on; plot(waves(2).x, waves(2).Dr)
xlim([0,500]); ylim([0 400]);
ylabel('D_r (W/m^2)')


subplot(5,1,5); plot(waves(1).x,waves(1).z,'k')
hold on; plot(waves(2).x, waves(2).z, 'k')
hold on; plot(waves(1).x,Zeta*ones(size(x)),'-.')
hold on; plot(waves(2).x,Zeta*ones(size(x)),'-.')
ylabel('zb (m)'); xlabel('x (m)')
xlim([0,500]);

