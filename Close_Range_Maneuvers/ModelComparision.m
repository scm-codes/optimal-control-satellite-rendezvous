% Sai Charan Malladi
% AE16B029 IIT Madras
% 2/1/2019

% This is the MAIN script for Model Comparision. The purpose of this code 
% is to evaluate several different methods of Model and implement control 
%
% Sources: 
%         -Howard D Curtis, Orbital Mechanics for Engineering Students
%         -Ryan Sherrill, Dynamics and Control of Satellite Relative Motion
%          in Elliptic Orbits using Lyapunov-Floquet Theory
%
% Credit to Quinn Kostelecky for the following function(s):
%     -PlanetaryConstants.m

% Here i have considered intial conditions to be given in only Classical
% Orbital Elements.

%% Begin 

clc; clear all;
%import the classical orbital elements of chaser/Deputy
CH = InitialClassicalOrbitalElements('ChaserCase1');
%import the classical orbital elments of Target/Chief
TR = InitialClassicalOrbitalElements('TargetCase1');

global mu;
mu = 3.986004415e5; 
m = 100;                                  % mass of the chief satellite
h_TR = sqrt(TR.a*mu*(1-TR.e^2));          % Angular momentum per unit mass of Target Satellite
T_TR = 2*pi*sqrt(TR.a^3/mu);              % Time Period of target Satellite
n_TR = 2*pi/T_TR;                         % Mean motion of target
h_CH = sqrt(CH.a*mu*(1-CH.e^2));          % Angular momentum of Chaser Satellite

%% Define Range of Times to Perform Analysis
numOrbits = 1;                              % number of orbits
t0 = 0;                                     % initial time               
tf = numOrbits*T_TR;                        % final time                  
numpoints = 20;                             % number of points to simulate 

% define time range
t_span = linspace(t0,tf,numpoints); 

%% Initial Conditions Determination 
% Initail Conditions(SV - statevector [x,y,z,xdot,ydot,zdot]; r - [x,y,z];
% v - [xdot, ydot, zdot]) are determined in the ECI frame from the given
% Classical Orbital Elements. 

%ECI fame
[r0_TR_ECI,v0_TR_ECI] = SV_from_COE(h_TR,TR.e,TR.raan,TR.i,TR.ap,TR.theta);
[r0_CH_ECI,v0_CH_ECI] = SV_from_COE(h_CH,CH.e,CH.raan,CH.i,CH.ap,CH.theta);

%Relative Position and Velocity of Chaser with respect to target in the ECI
%frame
rho0_ECI = r0_CH_ECI - r0_TR_ECI;
rhodot0_ECI = v0_CH_ECI - v0_TR_ECI;

%The relative positon and velocity are found in RSW frame.
[rho0_RSW, rhodot0_RSW] = RSW_from_SV(r0_TR_ECI,v0_TR_ECI,rho0_ECI,rhodot0_ECI);

%% Dynamic models in Relative frame attached to the target

%~~~~~~~~~~~~~~~~~~~~~ 1. Two-Body Propagation [TB]~~~~~~~~~~~~~~~~~~~~~~~~
%(Most Accuarate from which we derive Non linear equatuions of relative motion)
 
% Finding True Anomalies for every time step
% Estimating time since periapsis for the given initial conditions
tp_TR = TET_to_TP(TR.theta,TR.e,TR.a);
tp_CH = TET_to_TP(CH.theta,CH.e,CH.a);

%Define arrays with initial conditions
theta_TR = zeros(numpoints,1);
thetadot_TR = zeros(numpoints,1);
thetadotdot_TR = zeros(numpoints,1);
theta_CH = zeros(numpoints,1);
thetadot_CH = zeros(numpoints,1);
thetadotdot_CH = zeros(numpoints,1);


for i = 1:1:numpoints
   % Calculate true anomaly and Mean anamoly of target and chaser at every time step
   [theta_TR(i),M_TR(i)] = TP_to_TET(tp_TR+t_span(i),TR.e,TR.a);
   [theta_CH(i),M_CH(i)] = TP_to_TET(tp_CH+t_span(i),CH.e,CH.a);
end

for i = 1:1:numpoints
% Compute the derivaties of theta at each time step
[thetadot_TR(i),thetadotdot_TR(i)] = DERIVATIVES_of_THETA(theta_TR(i),TR.e,TR.a,TR.raan,TR.i,TR.ap,h_TR);
[thetadot_CH(i),thetadotdot_CH(i)] = DERIVATIVES_of_THETA(theta_CH(i),CH.e,CH.a,CH.raan,CH.i,CH.ap,h_CH);
end

%Define arrays
r_TR_TB_ECI = zeros(3,numpoints);
v_TR_TB_ECI = zeros(3,numpoints);
r_CH_TB_ECI = zeros(3,numpoints);1
v_CH_TB_ECI = zeros(3,numpoints);

rho_TB_ECI = zeros(3,numpoints);
rhodot_TB_ECI = zeros(3,numpoints);
rho_TB_RSW = zeros(3,numpoints);
rhodot_TB_RSW = zeros(3,numpoints);

for i = 1:1:numpoints
   % Estimate the state vector of target and chaser at each time step in
   % ECI frame from classical orbital elements
   [r_TR_TB_ECI(:,i),v_TR_TB_ECI(:,i)] = SV_from_COE(h_TR,TR.e,TR.raan,TR.i,TR.ap,theta_TR(i));
   [r_CH_TB_ECI(:,i),v_CH_TB_ECI(:,i)] = SV_from_COE(h_CH,CH.e,CH.raan,CH.i,CH.ap,theta_CH(i));
   
   % Find rho, rhodot in ECI frame
   rho_TB_ECI(:,i) = r_CH_TB_ECI(:,i) - r_TR_TB_ECI(:,i);
   rhodot_TB_ECI(:,i) = v_CH_TB_ECI(:,i) - v_TR_TB_ECI(:,i); 
   
   % Converting rho, rhodot from ECI to RSW
   [rho_TB_RSW(:,i),rhodot_TB_RSW(:,i)] = RSW_from_SV(r_TR_TB_ECI(:,i),v_TR_TB_ECI(:,i),rho_TB_ECI(:,i),rhodot_TB_ECI(:,i));
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~ 2. Linear Equtions of Relative Motion (LERM)~~~~~~~~~~~~~~~~
% LERM are derived by approximating that the distance between Chaser and
% Target is negligible when compared with distance of Target from ECI
% frame.

% Using ODE45 to solve for the time varying coefficients of linear
% equations.
opts = odeset('RelTol',1e-10,'AbsTol',1e-20);
[t,X] = ode45(@(t,X) LERM(t,X,tp_TR,TR.e,TR.a,TR.raan,TR.i,TR.ap,h_TR),[t0,tf],[rho0_RSW;rhodot0_RSW],opts);

%Define arrays
rho_LERM_RSW = zeros(numpoints,3);
rhodot_LERM_RSW = zeros(numpoints,3);

% interpolate for required number of points to compare
rho_LERM_RSW(:,:) = interp1(t,X(:,1:3),t_span(:));
rhodot_LERM_RSW(:,:) = interp1(t,X(:,4:6),t_span(:));

% Convert row to colounm vectors
rho_LERM_RSW = transpose(rho_LERM_RSW);
rhodot_LERM_RSW = transpose(rhodot_LERM_RSW);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~3. Hill - Clohessy - Wiltshire Equatios (HCW)~~~~~~~~~~~~~~~~~
% HCW equations are a special case of LERM where eccentricity e = 0. 
% Explicit solution for HCW exist. 

%Define array
rho_HCW_RSW = zeros(3,numpoints);
rhodot_HCW_RSW = zeros(3,numpoints);

%Using the HCW function to solve
for i = 1:1:numpoints
[rho_HCW_RSW(:,i),rhodot_HCW_RSW(:,i)] = HCW(rho0_RSW,rhodot0_RSW,n_TR,t_span(i)); 
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~4. Tschauner - Hempel Model (TH)~~~~~~~~~~~~~~~~~~~~~~
% TH equations are derived from LERM where the state matrix A and state
% vector x of x' = Ax are transformed to new state vector by coordinate
% scaling and the new state matrix A has independent variable as true
% anomaly(theta) unlike LERM where it is time elapsed.
% Refer Page No 20,21,22,23 in Sherrils disseration for more details

%Define array
rho_TH_RSW = zeros(3,numpoints);
rhodot_TH_RSW = zeros(3,numpoints);

% Use the TH fuction to solve
   [rho_TH_RSW,rhodot_TH_RSW] = TH(rho0_RSW,rhodot0_RSW,theta_TR,TR.a,TR.e,M_TR,h_TR,n_TR); 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~5. Virtual Cheif (VC)~~~~~~~~~~~~~~~~~~~~~~~~~~
% This is also called the dual deputy approach where the target (lets refer
% it as cheif here) and the chaser (lets refer it as deputy here) are
% considered moving with respect to a virtual cheif which has the same
% classical orbital elements as the real cheif with e = 0, so our HCW
% equations derived are more consistent with this frame. The actual
% relative position the difference between relative positions of actual
% cheif and depty using HCW equations wrt to virtual chief.
% For futher simplification of solution refer to pages 33, 34, 35 of
% sherril's disseration.

% Define arrays
rho_VC_RSW = zeros(3,numpoints);
rhodot_VC_RSW = zeros(3,numpoints);

% Use the function VC to solve
[rho_VC_RSW,rhodot_VC_RSW] = VC(theta_TR,M_TR,rho0_RSW,rhodot0_RSW,n_TR,thetadot_TR,t_span);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~6. Model in Felisiak's Paper on NMPC (FM)~~~~~~~~~~~~~~~~

% % Using ODE45 to solve for the time varying coefficients
% opts = odeset('RelTol',1e-10,'AbsTol',1e-20);
% [t,X] = ode45(@(t,X) FM(t,X,tp_TR,TR.e,TR.a,TR.raan,TR.i,TR.ap,h_TR),[t0,tf],[rho0_RSW;rhodot0_RSW],opts);
% 
% %Define arrays
% rho_FM_RSW = zeros(numpoints,3);
% rhodot_FM_RSW = zeros(numpoints,3);
% 
% % interpolate for required number of points to compare
% rho_FM_RSW(:,:) = interp1(t,X(:,1:3),t_span(:));
% rhodot_FM_RSW(:,:) = interp1(t,X(:,4:6),t_span(:));
% 
% % Convert row to colounm vectors
% rho_FM_RSW = transpose(rho_FM_RSW);
% rhodot_FM_RSW = transpose(rhodot_FM_RSW);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Evaluate Errors
% Find difference between approximates and "truth"
rhoerr_HCW_RSW = rho_HCW_RSW - rho_TB_RSW;
rhoerr_LERM_RSW = rho_LERM_RSW - rho_TB_RSW;
rhoerr_TH_RSW = rho_TH_RSW - rho_TB_RSW;
rhoerr_VC_RSW = rho_VC_RSW - rho_TB_RSW;
% rhoerr_FM_RSW = rho_FM_RSW - rho_TB_RSW;

% Find vector sum of errors
rhoerr_LERM_mag = zeros(1,numpoints);
rhoerr_HCW_mag = zeros(1,numpoints);
rhoerr_TH_mag = zeros(1,numpoints);
rhoerr_VC_mag = zeros(1,numpoints);
% rhoerr_FM_mag = zeros(1,numpoints);

for i = 1:numpoints
    rhoerr_LERM_mag(i) = norm(rhoerr_LERM_RSW(:,i));
    rhoerr_HCW_mag(i) = norm(rhoerr_HCW_RSW(:,i));
    rhoerr_TH_mag(i) = norm(rhoerr_TH_RSW(:,i));
    rhoerr_VC_mag(i) = norm(rhoerr_VC_RSW(:,i));
%     rhoerr_FM_mag(i) = norm(rhoerr_FM_RSW(:,i));
end

% Find RMS error
rhoerr_LERM = rms(rhoerr_LERM_mag);
rhoerr_HCW = rms(rhoerr_HCW_mag);
rhoerr_TH = rms(rhoerr_TH_mag);
rhoerr_VC = rms(rhoerr_VC_mag);
% rhoerr_FM = rms(rhoerr_FM_mag);

fprintf('HCW RMS: %f\n',rhoerr_HCW)
fprintf('TH RMS: %f\n',rhoerr_TH)
fprintf('LERM RMS: %f\n',rhoerr_LERM)
fprintf('VC RMS: %f\n',rhoerr_VC)
% fprintf('FM RMS: %f\n',rhoerr_FM)

%% Plot
set(groot,'DefaultAxesColorOrder',[0 0 1; 1 0 0; 0 0 0],'DefaultLineLineWidth',2);
set(groot,'DefaultAxesFontSize',20,'defaultAxesFontName', 'Times','defaultTextFontName', 'Times')

figure;
hold on
grid on
view(3);
p1 = plot3(r_TR_TB_ECI(1,:),r_TR_TB_ECI(2,:),r_TR_TB_ECI(3,:),'ro-','LineWidth',2);
p2 = plot3(r_CH_TB_ECI(1,:),r_CH_TB_ECI(2,:),r_CH_TB_ECI(3,:),'kx-','LineWidth',2);
xlabel('i (in Km)');
ylabel('j (in Km)');
zlabel('k (in Km)');
title('Trajectories of Chief and Deputy in ECI Frame');
legend([p1,p2],'Chief','Deputy')
hold off

figure;
hold on 
grid on
p1 = plot(t_span(:),rhoerr_LERM_mag(:),'bo-','LineWidth',2);
p2 = plot(t_span(:),rhoerr_HCW_mag(:),'mx-','LineWidth',2);
p3 = plot(t_span(:),rhoerr_TH_mag(:),'k+-','LineWidth',2);
p4 = plot(t_span(:),rhoerr_VC_mag(:),'r*-','LineWidth',2);
% p5 = plot(t_span(:),rhoerr_FM_mag(:),'k-','LineWidth',2);
xlabel('Time (in s)')
ylabel('Relative Position error at each time step (in Km)');
legend([p1,p2,p3,p4],'LERM','HCW','TH','VC');
title('Comparision of error at each time step for one Cheif Time Period');
hold off

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~