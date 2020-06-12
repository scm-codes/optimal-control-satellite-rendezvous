% Sai Charan Malladi 
% AE16B029 IIT Madras 
% BTP Work
% 2/4/2019
%
% This code uses the TH model and uses Lyapunov-Floquet theory to change
% the periodic system into invariant system and then uses lqr for control.
%
% Here i have considered intial conditions to be given in only Classical
% Orbital Elements.

clc; clear all;
%import the classical orbital elements of chaser
CH = InitialClassicalOrbitalElements('ChaserCase1');
%import the classical orbital elments of Target
TR = InitialClassicalOrbitalElements('TargetCase1');

global mu;
mu = 3.986004415e5; 
m = 100; % mass of the chief satellite
h_TR = sqrt(TR.a*mu*(1-TR.e^2)); % Angular momentum per unit mass of Target Satellite
T_TR = 2*pi*sqrt(TR.a^3/mu); % Time Period of target Satellite
n = 2*pi/T_TR; % Mean motion of target
h_CH = sqrt(CH.a*mu*(1-CH.e^2));% Angular momentum of Chaser Satellite
p_TR = h_TR^2/mu;
% c = sqrt(1+3*EARTH.J2*EARTH.Re^2/(8*TR.a^2)*(1+3*cos(2*TR.i))); % A constant entity used in modelling J2

%% Define Range of Times to Perform Analysis
numOrbits = 1;                           % number of orbits
t0 = 0;                                     % initial time               
tf = numOrbits*T_TR;                        % final time                   
Ts = 1;                                     % sampling

% define time range
t_span = t0:Ts:tf;
numpoints = length(t_span);

%% Initial Conditions Determination
% Initail Conditions(SV - statevector [x,y,z,xdot,ydot,zdot]; r - [x,y,z]; v - [xdot, ydot, zdot]) are
% determined in the ECI frame from the given Classical Orbital Elements.

%ECI fame
[r0_TR_ECI,v0_TR_ECI] = SV_from_COE(h_TR,TR.e,TR.raan,TR.i,TR.ap,TR.theta);
[r0_CH_ECI,v0_CH_ECI] = SV_from_COE(h_CH,CH.e,CH.raan,CH.i,CH.ap,CH.theta);

%Relative Position and Velocity of Chaser with respect to target in the ECI
%frame
rho0_ECI = r0_CH_ECI - r0_TR_ECI;
rhodot0_ECI = v0_CH_ECI - v0_TR_ECI;

%The relative positon and velocity are found in RSW frame.
[rho0_RSW, rhodot0_RSW] = RSW_from_SV(r0_TR_ECI,v0_TR_ECI,rho0_ECI,rhodot0_ECI);

%% TH

A =  [0            0      0            1       0       0
      0            0      0            0       1       0
      0            0      0            0       0       1
      3*n^2        0      0            0       2*n     0
      0            0      0           -2*n     0       0
      0            0     -1*n^2        0       0       0];
  
% Control Matrix using 
B =[0 0 0
    0 0 0
    0 0 0 
    1 0 0
    0 1 0
    0 0 1];

% Q = [1 0 0 0 0 0
%     0 1 0 0 0 0
%     0 0 1 0 0 0
%     0 0 0 1/n^2 0 0
%     0 0 0 0 1/n^2 0
%     0 0 0 0 0 1/n^2];
% 
% R = 100*[1/n^4 0 0 
%     0 1/n^4 0
%     0 0 1/n^4];

Q = [1 0 0 0 0 0
    0 1 0 0 0 0
    0 0 1 0 0 0
    0 0 0 1 0 0
    0 0 0 0 1 0
    0 0 0 0 0 1];

R = 10^4*[1 0 0 
    0 1 0
    0 0 1];

[K,S,E] = lqr(A,B,Q,R);

%% For computation of A(t)

% Finding True Anomalies for every time step
% Estimating time since periapsis for the given initial conditions
tp_TR = TET_to_TP(TR.theta,TR.e,TR.a);
tp_CH = TET_to_TP(CH.theta,CH.e,CH.a);
[theta0,M0] = TP_to_TET(tp_TR,TR.e,TR.a);

% P0 chosen for periapse matching transformation
e = TR.e;
P0 = [2*(1-e^2)^(5/2)/((e+1)^3*(e+2)) 0 0 0 ((1-e^2)^(5/2))/(n*(e+1)^3*(e+2)) - 1/(2*n) 0
    0 1 0 0 0 0
    0 0 1 0 0 0
    0 e*h_TR*(e+1)/p_TR^2 0 n*p_TR^2*(e+1)^2/(h_TR*(1-e^2)^(5/2)) 0 0
    0 0 0 0 h_TR*(e+1)*(e+2)/(2*n*p_TR^2) 0
    0 0 0 0 0 1];

%  
K0 = M0/(1-e^2)^(3/2);
PSI0 = getPSI(K0,theta0,e);
T0 = [(1+e*cos(theta0))*eye(3),             zeros(3);
    -e*sin(theta0)*eye(3),      p_TR^2/(h_TR*(1+e*cos(theta0)))*eye(3)];

% U
U = zeros(3,numpoints);

% simulation
X(:,1) = [rho0_RSW;rhodot0_RSW];

for i = 1:1:numpoints
%Compute theta and its derivative at each time as the coefficients of the
%linear equations vary with time.
t = t_span(i);
[theta,M] = TP_to_TET(t.*1+tp_TR,TR.e,TR.a);
r = h_TR^2/mu/(1+e*cos(theta));
[thetadot,thetadotdot] = DERIVATIVES_of_THETA(theta,TR.e,TR.a,TR.raan,TR.i,TR.ap,h_TR);

    a41 = thetadot^2 + 2*mu/r^3;
    a42 = thetadotdot;
    a45 = 2*thetadot;
    a51 = -thetadotdot;
    a52 = thetadot^2 - mu/r^3;
    a54 = - 2*thetadot;
    a63 = - mu/r^3;
    
    %%%% Constructing LERM Model
     A = [0 0 0 1 0 0
    0 0 0 0 1 0
    0 0 0 0 0 1
    a41 a42 0 0 a45 0 
    a51 a52 0 a54 0 0
    0 0 a63 0 0 0];
    
    % Method described by Sherill in his paper
    invphiHCW = expm(-A.*t);

    Z = M/(1-e^2)^(3/2);
    
    PSI = getPSI(Z,theta,e);

    % find T
    T = [(1+e*cos(theta))*eye(3),            zeros(3);
    -e*sin(theta)*eye(3),       p_TR^2/(h_TR*(1+e*cos(theta)))*eye(3)];

    % find state at new time
    phiLERM = inv(T)*PSI*inv(PSI0)*T0;
    
    % Find P at new time
    P = phiLERM*P0*invphiHCW;
    K_bar = B'*P*B*K*inv(P);
    dX = (A - B*K_bar)*X(:,i);
    U(:,i) = -K_bar*X(:,i);
    X(:,i+1) = X(:,i) + dX;
    eig(A - B*K_bar)
    
    % Stop Condition
    if rms(X(1:3,i)) <= 0.0005
       endtime = i;
       break;
    else
       endtime = i;
    end  
end


%% Plots
set(groot,'DefaultAxesColorOrder',[0 0 1; 1 0 0; 0 0 0],'DefaultLineLineWidth',2);
set(groot,'DefaultAxesFontSize',20,'defaultAxesFontName', 'Times','defaultTextFontName', 'Times')

figure;
hold on
grid on
view(3);
p0 = plot3(X(1,1),X(2,1),X(3,1),'mo','LineWidth',2);
p1 = plot3(X(1,1:endtime),X(2,1:endtime),X(3,1:endtime),'b-.','LineWidth',2);
title('Relative Motion Trajectory of Satellites');
xlabel('x (in Km)');
ylabel('y (in Km)');
zlabel('z (in Km)');
legend('Start','LF')
hold off

figure;
suptitle('Relative Position');
subplot(3,1,1)
hold on
grid on
p1 = plot(t_span(1:endtime),X(1,1:endtime),'b-.','LineWidth',2);
legend('LF')
xlabel('t (in s)');
ylabel('x (in Km)');
hold off

subplot(3,1,2)
hold on
grid on
p1 = plot(t_span(1:endtime),X(2,1:endtime),'b-.','LineWidth',2);
xlabel('t (in s)');
ylabel('y (in Km)');
hold off

subplot(3,1,3)
hold on
grid on
p1 = plot(t_span(1:endtime),X(3,1:endtime),'b-.','LineWidth',2);
xlabel('t (in s)');
ylabel('z (in Km)');
hold off

figure;
suptitle('Control Input');
subplot(3,1,1)
hold on
grid on
p1 = plot(t_span(1:endtime),U(1,1:endtime),'b-.','LineWidth',2);
legend('LF')
xlabel('t (in s)');
ylabel('Ux (in KN)');
hold off

subplot(3,1,2)
hold on
grid on
p1 = plot(t_span(1:endtime),U(2,1:endtime),'b-.','LineWidth',2);
xlabel('t (in s)');
ylabel('Uy (in KN)');
hold off

subplot(3,1,3)
hold on
grid on
p1 = plot(t_span(1:endtime),U(3,1:endtime),'b-.','LineWidth',2);
xlabel('t (in s)');
ylabel('Uz (in KN)');
hold off