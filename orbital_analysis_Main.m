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

% Assumptions:
% Sun is approximated as a point mass
% Other planets are not exhibiting torques on the station

% 3 Different shapes for stations and accompannying moments of inertia
% Cylinder, Ring, and Growth Adaptable Artificial Gravity Space Habitat

% Mass values for small, medium, and large habitats
m_s = 1e5; m_m = 1e9; m_l = 1e14;               % metric tons
M = (1e3)*[m_s m_m m_l];

% Cylinder dimensions for small, medium, and large habitats
R_s = 100; R_m = 4e3; R_l = 3.2e6;              % m
Rc = [R_s R_m R_l];
L_s = 500; L_m = 2e4; L_l = 1.6e7;              % m
Lc = [L_s L_m L_l];

% Ring dimensions for small, medium, and large habitats
R_1s = 280; R_1m = 1.12e4; R_1l = 9.95e6;       % m
Rr_1 = [R_1s R_1m R_1l];
R_2s = 300; R_2m = 1.2e4; R_2l = 1e7;           % m
Rr_2 = [R_2s R_2m R_2l];
h_s = 30; h_m = 1.2e3; h_l = 3.18e5;            % m
hr = [h_s h_m h_l];

% GAAGSH dimensions for small, medium, and large habitats
Mcyl = (3/5)*M; Mcone = (1/5)*M;                % kg
R_s = 224; R_m = 8.96e3; R_l = 7.168e6;         % m
Rsh = [R_s R_m R_l];
Lt_s = 259; Lt_m = 1.036e4; Lt_l = 8.288e6;     % m
Lt = [Lt_s Lt_m Lt_l];
Lh_s = 517; Lh_m = 2.068e4; Lh_l = 1.6544e7;    % m
Lh = [Lh_s Lh_m Lh_l];

% Cylinder Principal Moments of Inertia 
% Note: Each row is a different size set of moment of inerties, while the
% columns are the different prinicpal moments of inertia per axes
Ic = [];
for i = 1:length(M)
    for j = 1:3
        if j ~= 3
            Ic(i,j) = (1/12)*M(i)*Lc(i)^2;
        else
            Ic(i,j) = M(i)*Rc(i)^2;
        end
    end
end

% Ring Principal Moments of Intertia
Ir = [];
for i = 1:length(M)
    for j = 1:3
        if j ~= 3
            Ir(i,j) = (1/12)*M(i)*(3*(Rr_1(i)^2 + Rr_2(i)^2) + hr(i)^2); 
        else
            Ir(i,j) = (1/2)*M(i)*(Rr_1(i)^2+Rr_2(i)^2);
        end
    end
end

% Space Habitat Principal Moments of Intertia
Ish = [];
for i = 1:length(M)
    for j = 1:3
        if j ~= 3
            Ish(i,j) = (1/12)*Mcyl(i)*Lh(i)^2 + (1/2)*Mcone(i)*(Rsh(i)^2 + (1/2)*(Lh(i)^2 + Lt(i)^2 + Lt(i)) + Lh(i) - Lh(i)*Lt(i)); 
        else
            Ish(i,j) = (Mcyl(i)+Mcone(i))*(3/2)*Rsh(i)^2;
        end
    end
end

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

%% Cylinder Analysis (Gravity Gradient only)

% Setting initial angles and transforming them to quaternions
ang0 = [0 pi/2 0];
R0 = DCM321(ang0);
q0 = r_q(R0);

% Setting matrix with initial angular rates (rows are different size
% habitats, columns are the body axis rates)
w0 = [-n*cos(ang0(2)) 0 n*sin(ang0(2))];
s0 = [q0 w0(1,:)];

k = 2;
tic
fprintf('Cylinder rotational analysis started... \n')
qc = []; wc = [];
for i = 1:length(M)
    tic
    [t,s] = ode45(@(t,s) DiffEQ_Q(t,s,Ic(i,1:3)),tspan,s0,1e-9);
    toc
    qc{i} = s(:,1:4); wc{i} = s(:,5:7);
    psi_c = zeros(length(M),length(qc)); theta_c = zeros(length(M),length(qc)); phi_c = zeros(length(M),length(qc));
    for j = 1:length(qc{i})
        qc{i}(j,:) = qc{i}(j,:)/norm(qc{i}(j,:));
        [psi_c(i,j) theta_c(i,j) phi_c(i,j)] = Qto321Eul(qc{i}(j,:));
    end
    
    figure(k)
    subplot(3,1,1)
    plot(t/T_y,rad2deg(psi_c(i,:)))
    title('\Psi')
    ylabel(['Yaw (' char(176) ')'])
    subplot(3,1,2)
    plot(t/T_y,rad2deg(theta_c(i,:)))
    title('\Theta')
    ylabel(['Pitch (' char(176) ')'])
    subplot(3,1,3)
    plot(t/T_y,rad2deg(phi_c(i,:)))
    title('\Phi')
    ylabel(['Roll (' char(176) ')'])
    xlabel('Time (Years)')
    if i == 1
        fprintf('Small structure analysis complete.\n')
    elseif i == 2
        fprintf('Medium structure analysis complete.\n')
    elseif i == 3
        fprintf('Large structure analysis complete.\n')
    end
    k = k+1;
end
toc
fprintf('Cylinder rotational analysis complete!\n\n')

%% Ring Analysis (Gravity Gradient only)

% Setting initial angles and transforming them to quaternions
ang0 = [0 -pi/4 0];
R0 = DCM321(ang0);
q0 = r_q(R0);

% Setting matrix with initial angular rates (rows are different size
% habitats, columns are the body axis rates)
w0 = [-n*cos(ang0(2)) 0 n*sin(ang0(2))];
s0 = [q0 w0];

tic
fprintf('Ring rotational analysis started... \n')
qr = []; wr = [];
for i = 1:length(M)
    tic
    [t,s] = ode45(@(t,s) DiffEQ_Q(t,s,Ir(i,1:3)),tspan,s0,1e-9);
    toc
    qr{i} = s(:,1:4); wr{i} = s(:,5:7);
    psi_r = zeros(length(M),length(qr)); theta_r = zeros(length(M),length(qr)); phi_r = zeros(length(M),length(qr));
    for j = 1:length(qr{i})
        qr{i}(j,:) = qr{i}(j,:)/norm(qr{i}(j,:));
        [psi_r(i,j) theta_r(i,j) phi_r(i,j)] = Qto321Eul(qr{i}(j,:));
    end
    
    figure(k)
    subplot(3,1,1)
    plot(t/T_y,rad2deg(psi_r(i,:)))
    title('\Psi')
    ylabel(['Yaw (' char(176) ')'])
    subplot(3,1,2)
    plot(t/T_y,rad2deg(theta_r(i,:)))
    title('\Theta')
    ylabel(['Pitch (' char(176) ')'])
    subplot(3,1,3)
    plot(t/T_y,rad2deg(phi_r(i,:)))
    title('\Phi')
    ylabel(['Roll (' char(176) ')'])
    xlabel('Time (Years)')
    if i == 1
        fprintf('Small structure analysis complete.\n')
    elseif i == 2
        fprintf('Medium structure analysis complete.\n')
    elseif i == 3
        fprintf('Large structure analysis complete.\n')
    end
    k = k+1;
end
toc
fprintf('Ring rotational analysis complete!\n\n')

%% Space Habitat Analysis (Gravity Gradient only)

% Setting initial angles and transforming them to quaternions
ang0 = [0 0 0];
R0 = DCM321(ang0);
q0 = r_q(R0);

% Setting matrix with initial angular rates (rows are different size
% habitats, columns are the body axis rates)
w0 = [-n*cos(ang0(2)) 0 n*sin(ang0(2))];
s0 = [q0 w0];

tic
fprintf('Space Habitat rotational analysis started... \n')
qs = []; ws = [];
for i = 1:length(M)
    tic
    [t,s] = ode45(@(t,s) DiffEQ_Q(t,s,Ish(i,1:3)),tspan,s0,1e-9);
    toc
    qs{i} = s(:,1:4); ws{i} = s(:,5:7);
    psi_s = zeros(length(M),length(qs)); theta_s = zeros(length(M),length(qs)); phi_s = zeros(length(M),length(qs));
    for j = 1:length(qs{i})
        qs{i}(j,:) = qs{i}(j,:)/norm(qs{i}(j,:));
        [psi_s(i,j) theta_s(i,j) phi_s(i,j)] = Qto321Eul(qs{i}(j,:));
    end
    
    figure(k)
    subplot(3,1,1)
    plot(t/T_y,rad2deg(psi_s(i,:)))
    title('\Psi')
    ylabel(['Yaw (' char(176) ')'])
    subplot(3,1,2)
    plot(t/T_y,rad2deg(theta_s(i,:)))
    title('\Theta')
    ylabel(['Pitch (' char(176) ')'])
    subplot(3,1,3)
    plot(t/T_y,rad2deg(phi_s(i,:)))
    title('\Phi')
    ylabel(['Roll (' char(176) ')'])
    xlabel('Time (Years)')
    if i == 1
        fprintf('Small structure analysis complete.\n')
    elseif i == 2
        fprintf('Medium structure analysis complete.\n')
    elseif i == 3
        fprintf('Large structure analysis complete.\n')
    end
    k = k+1;
end
toc
fprintf('Space Habitat rotational analysis complete!\n\n')

%% Cylinder Analysis (Gravity Gradient and Artificial Gravity Rotation)

% Setting initial angles and transforming them to quaternions
ang0 = [0 pi/2 0];
R0 = DCM321(ang0);
q0 = r_q(R0);

% Setting matrix with initial angular rates (rows are different size
% habitats, columns are the body axis rates)
w0 = [0.221-n*cos(ang0(2)) 0 n*sin(ang0(2));...
        (4.95e-2)-n*cos(ang0(2)) 0 n*sin(ang0(2));...
        (1.75e-3)-n*cos(ang0(2)) 0 n*sin(ang0(2))];
% Setting initial conditions for each case
s0 = [q0 w0(1,:);...
        q0 w0(2,:);...
        q0 w0(3,:)];

tic
fprintf('Cylinder rotational analysis started... \n')
qc2 = []; wc2 = [];
for i = 1:length(M)
    tic
    [t,s] = ode45(@(t,s) DiffEQ_Q(t,s,Ic(i,1:3)),tspan,s0(i,:),1e-9);
    toc
    qc2{i} = s(:,1:4); wc2{i} = s(:,5:7);
    psi_c2 = zeros(length(M),length(qc2)); theta_c2 = zeros(length(M),length(qc2)); phi_c2 = zeros(length(M),length(qc2));
    for j = 1:length(qc2{i})
        qc2{i}(j,:) = qc2{i}(j,:)/norm(qc2{i}(j,:));
        [psi_c2(i,j) theta_c2(i,j) phi_c2(i,j)] = Qto321Eul(qc2{i}(j,:));
    end
    
    figure(k)
    subplot(3,1,1)
    plot(t/T_y,rad2deg(psi_c2(i,:)))
    title('\Psi')
    ylabel(['Yaw (' char(176) ')'])
    subplot(3,1,2)
    plot(t/T_y,rad2deg(theta_c2(i,:)))
    title('\Theta')
    ylabel(['Pitch (' char(176) ')'])
    subplot(3,1,3)
    plot(t/T_y,rad2deg(phi_c2(i,:)))
    title('\Phi')
    ylabel(['Roll (' char(176) ')'])
    xlabel('Time (Years)')
    if i == 1
        fprintf('Small structure analysis complete.\n')
    elseif i == 2
        fprintf('Medium structure analysis complete.\n')
    elseif i == 3
        fprintf('Large structure analysis complete.\n')
    end
    k = k+1;
end
toc
fprintf('Cylinder rotational analysis with artificial gravity complete!\n\n')

%% Ring Analysis (Gravity Gradient and Artificial Gravity Rotation)

% Setting initial angles and transforming them to quaternions
ang0 = [0 -pi/4 0];
R0 = DCM321(ang0);
q0 = r_q(R0);

% Setting matrix with initial angular rates (rows are different size
% habitats, columns are the body axis rates)
w0 = [0.187-n*cos(ang0(2)) 0 n*sin(ang0(2));...
        (2.96e-2)-n*cos(ang0(2)) 0 n*sin(ang0(2));...
        (9.92e-4)-n*cos(ang0(2)) 0 n*sin(ang0(2))];

% Setting initial conditions for each case
s0 = [q0 w0(1,:);...
        q0 w0(2,:);...
        q0 w0(3,:)];

tic
fprintf('Ring rotational analysis started... \n')
qr2 = []; wr2 = [];
for i = 1:length(M)
    tic
    [t,s] = ode45(@(t,s) DiffEQ_Q(t,s,Ir(i,1:3)),tspan,s0(i,:),1e-9);
    toc
    qr2{i} = s(:,1:4); wr2{i} = s(:,5:7);
    psi_r2 = zeros(length(M),length(qr2)); theta_r2 = zeros(length(M),length(qr2)); phi_r2 = zeros(length(M),length(qr2));
    for j = 1:length(qr2{i})
        qr2{i}(j,:) = qr2{i}(j,:)/norm(qr2{i}(j,:));
        [psi_r2(i,j) theta_r2(i,j) phi_r2(i,j)] = Qto321Eul(qr2{i}(j,:));
    end
    
    figure(k)
    subplot(3,1,1)
    plot(t/T_y,rad2deg(psi_r2(i,:)))
    title('\Psi')
    ylabel(['Yaw (' char(176) ')'])
    subplot(3,1,2)
    plot(t/T_y,rad2deg(theta_r2(i,:)))
    title('\Theta')
    ylabel(['Pitch (' char(176) ')'])
    subplot(3,1,3)
    plot(t/T_y,rad2deg(phi_r2(i,:)))
    title('\Phi')
    ylabel(['Roll (' char(176) ')'])
    xlabel('Time (Years)')
    if i == 1
        fprintf('Small structure analysis complete.\n')
    elseif i == 2
        fprintf('Medium structure analysis complete.\n')
    elseif i == 3
        fprintf('Large structure analysis complete.\n')
    end
    k = k+1;
end
toc
fprintf('Ring rotational analysis with artificial gravity complete!\n\n')

%% Space Habitat Analysis (Gravity Gradient and Artificial Gravity Rotation)

% Setting initial angles and transforming them to quaternions
ang0 = [0 0 0];
R0 = DCM321(ang0);
q0 = r_q(R0);

% Setting matrix with initial angular rates (rows are different size
% habitats, columns are the body axis rates)
w0 = [0.221-n*cos(ang0(2)) 0 n*sin(ang0(2));...
        (3.50e-2)-n*cos(ang0(2)) 0 n*sin(ang0(2));...
        (1.24e-3)-n*cos(ang0(2)) 0 n*sin(ang0(2))];
% Setting initial conditions for each case
s0 = [q0 w0(1,:);...
        q0 w0(2,:);...
        q0 w0(3,:)];

tic
fprintf('Space Habitat rotational analysis started... \n')
qs2 = []; ws2 = [];
for i = 1:length(M)
    tic
    [t,s] = ode45(@(t,s) DiffEQ_Q(t,s,Ish(i,1:3)),tspan,s0(i,:),1e-9);
    toc
    qs2{i} = s(:,1:4); ws2{i} = s(:,5:7);
    psi_s2 = zeros(length(M),length(qs2)); theta_s2 = zeros(length(M),length(qs2)); phi_s = zeros(length(M),length(qs2));
    for j = 1:length(qs2{i})
        qs{i}(j,:) = qs{i}(j,:)/norm(qs{i}(j,:));
        [psi_s2(i,j) theta_s2(i,j) phi_s2(i,j)] = Qto321Eul(qs2{i}(j,:));
    end
    
    figure(k)
    subplot(3,1,1)
    plot(t/T_y,rad2deg(psi_s2(i,:)))
    title('\Psi')
    ylabel(['Yaw (' char(176) ')'])
    subplot(3,1,2)
    plot(t/T_y,rad2deg(theta_s2(i,:)))
    title('\Theta')
    ylabel(['Pitch (' char(176) ')'])
    subplot(3,1,3)
    plot(t/T_y,rad2deg(phi_s2(i,:)))
    title('\Phi')
    ylabel(['Roll (' char(176) ')'])
    xlabel('Time (Years)')
    if i == 1
        fprintf('Small structure analysis complete.\n')
    elseif i == 2
        fprintf('Medium structure analysis complete.\n')
    elseif i == 3
        fprintf('Large structure analysis complete.\n')
    end
    k = k+1;
end
toc
fprintf('Space Habitat rotational analysis with artificial gravity complete!\n\n')

%% Functions

function ds = DiffEQ(t,s)
    mass_S = 1.989e30; G = 6.67408e-11; mu_s = G*mass_S*1e-9;
    rv = s(1:3); vv = s(4:6);
    r = norm(rv);
    
    ds(1:3) = vv;
    ds(4:6) = -(mu_s/(r^3))*(rv);
    ds = ds';
end

function ds = DiffEQ2(t,s,I);
    mass_S = 1.989e30; G = 6.67408e-11; mu_s = G*mass_S*1e-9; AU = 1.495978707e8;
    rv = s(1:3); vv = s(4:6); angv = s(7:9); wv = s(10:12);
    phi = angv(1); theta = angv(2); psi = angv(3);
    It = I(1); Iz = I(3);
    r = norm(rv);
    Omega_sq = mu_s/(r)^3;
    c1 = cos(theta)*cos(psi);
    c2 = sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi);
    c3 = cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);
    Mx = 3*Omega_sq*(Iz-It)*c2*c3;
    My = 3*Omega_sq*(It-Iz)*c1*c3;
    
    ds(1:3) = vv;
    ds(4:6) = -(mu_s/(r^3))*(rv);
    ds(7) = wv(2)*sin(phi)/cos(theta) + wv(3) - sqrt(Omega_sq);
    ds(8) = wv(2)*cos(phi) - wv(3)*cos(phi);
    ds(9) = wv(1) + wv(2)*sin(phi)*tan(theta) + wv(3)*cos(phi)*tan(theta);
    ds(10) = 0;
    ds(11) = (1/It)*(My - (Iz-It)*wv(1)*wv(3));
    ds(12) = (1/It)*(Mx - (It-Iz)*wv(1)*wv(2));
    ds = ds';
    t
end

function ds = DiffEQ_Q(t,s,I);
    mass_S = 1.989e30; G = 6.67408e-11; mu_s = G*mass_S*1e-9; AU = 1.495978707e8;
    q = s(1:4); % q = q/norm(q);
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    wv = s(5:7);
    It = I(1); Iz = I(3);
    Omega_sq = mu_s/(AU)^3;
    c1 = q3^2 + q0^2 - q1^2 - q2^2;
    c2 = 2*(q0*q1-q3*q2);
    c3 = 2*(q0*q2+q3*q1);
    Mx = 3*Omega_sq*(Iz-It)*c2*c3;
    My = 3*Omega_sq*(It-Iz)*c1*c3;
    q_Mat = QuatMat(q);
    
    
    ds(1:4) = 0.5 * q_Mat * [0; wv];
    ds(5) = 0;
    ds(6) = (1/It)*(My - (Iz-It)*wv(1)*wv(3));
    ds(7) = (1/It)*(Mx - (It-Iz)*wv(1)*wv(2));
    ds = ds';
%     t
end


