clear, clc, close all

%% AERO 624 Celestial Mechanics - William Young
% 4/21/20
% Class Project
% Project Description: Analysis of large space habitats and their
% associated garvity gradient torques in the deep space (0.9 - 1.1 AU)
% regime

% Constants
mass_S = 1.989e30;                              % kg
G = 6.67408e-11;                                % m^3 / [kg * s^2]
mu_s = G*mass_S*1e-9;                           % km^3 * s^-2
AU = 1.495978707e8;                             % km
T_y = 365.25*3600*24;                           % s

% Space Station Orbital Regime
Orb_regime = [0.9 1.1];                         % AU

rv0 = [AU 0 0];
vv0 = [0 sqrt(mu_s/norm(rv0)) 0];
n = sqrt(mu_s/(AU)^3);

%% Trajectory Analysis

% Initial State
s0 = [rv0 vv0];

% Setting time spans
T = 2*pi/n;
tspan = linspace(0,T,1e6);
tspan1 = linspace(0,2*T,1e6);
tspan2 = [0 T]; tspan3 = [0 2*T]; tspan5 = [0 5*T];

% Attempting to normalize
tic
mu_sn = mus_s/AU;
r_vn = [1 0 0];
v_vn = [0 sqrt((mu_sn)/norm(r_vn)) 0];
n = sqrt(mu_sn/(AU)^3);
toc

tic
fprintf('Trajectory analysis started... \n')
[t,s] = ode45(@(t,s) DiffEQ(t,s),linspace(0,T,1e6),s0,1e-9);
r = s(:,1:3); r_sc = r/AU; v = s(:,4:6);
fprintf('Error with 1e6 timesteps for 1 period (1 year): %1.4e \n',norm(r(end,:)-r(1,:)))
[t,s] = ode45(@(t,s) DiffEQ(t,s),tspan1,s0,1e-9);
fprintf('Error with 1e6 timesteps for 2 periods (2 years): %1.4e \n',norm(s(round(length(s)/2),1:3)-s(1,1:3)))
[t,s] = ode45(@(t,s) DiffEQ(t,s),tspan2,s0,1e-9);
fprintf('Error with matlab set timesteps for 1 periods (1 year): %1.4e \n',norm(s(end,1:3)-s(1,1:3)))
[t,s] = ode45(@(t,s) DiffEQ(t,s),tspan3,s0,1e-9);
fprintf('Error with matlab set timesteps for 2 periods (2 years): %1.4e \n',norm(s(round(length(s)/2),1:3)-s(1,1:3)))
[t,s1] = ode45(@(t,s) DiffEQ(t,s),tspan5,s0,1e-9);
fprintf('Absolute error for 1 year: %4.4f %% \n',norm(r(end,:)-r(1,:))/AU*100)
[t,s] = ode45(@(t,s) DiffEQ(t,s),tspan1,s0,1e-9);
[t,s1] = ode45(@(t,s) DiffEQ(t,s),linspace(0,5*T,1e6),s0,1e-9);
fprintf('Absolute error for 2 years: %4.4f %% \n',norm(s(round(length(s)/2),1:3)-s(1,1:3))/AU*100)
fprintf('Absolute error for 5 years: %4.4f %% \n',norm(s1(round(length(s1)/5),1:3)-s1(1,1:3))/AU*100)
toc
fprintf('Trajectory Analysis complete! \n\n')

figure()
figure(1)
plot(s1(:,1)/AU,s1(:,2)/AU,s(:,1)/AU,s(:,2)/AU,r_sc(:,1),r_sc(:,2))
title('Orbit Plot')
legend('1 AU orbit propagated for 5 years','1 AU orbit propagated for 2 years','1 AU orbit propagated for 1 year')
xlabel('x (AU)')
ylabel('y (AU)')
zlabel('z (AU)')

%% Functions

function ds = DiffEQ(t,s)
    mass_S = 1.989e30; G = 6.67408e-11; mu_s = G*mass_S*1e-9;
    rv = s(1:3); vv = s(4:6);
    r = norm(rv);
    
    ds(1:3) = vv;
    ds(4:6) = -(mu_s/(r^3))*(rv);
    ds = ds';
end