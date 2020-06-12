% Sai Charan Malladi 
% AE16B029 IIT Madras 
% Winter Work
% 1/1/2019

% This code uses the HCW model and uses LQR
% approach to arrie at optimal gain and thus generating optimal
% trajectories in LVLH(RSW) frame

%Here i have considered intial conditions to be given in only Classical
%Orbital Elements.

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
h_CH = sqrt(CH.a*mu*(1-CH.e^2)); % Angular momentum of Chaser Satellite
% c = sqrt(1+3*EARTH.J2*EARTH.Re^2/(8*TR.a^2)*(1+3*cos(2*TR.i))); % A constant entity used in modelling J2

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

%% Define Range of Times to Perform Analysis
numOrbits = 1/4;                            % number of orbits
t0 = 0;                                     % initial time               
tf = numOrbits*T_TR;                        % final time                   
Ts = 1;                                     % sampling

% define time range
t_span = t0:Ts:tf;
numpoints = length(t_span);

%% Control Part Using LQR

%%%% Trajectory in LVLH frame with LQR control

% Initial HCW state space model
A = [0            0      0            1       0       0
      0            0      0            0       1       0
      0            0      0            0       0       1
     -3*n^2        0      0            0       2*n     0
      0            0      0           -2*n     0       0
      0            0     -1*n^2        0       0       0];
  
% Control Matrix using 
B = 1/m*[ 0 0 0
          0 0 0
          0 0 0 
          1 0 0
          0 1 0
          0 0 1];
C = eye(6);
D = 0.*eye(6,3);
sys = ss(A,B,C,D);

% If rank of the controllabity matrix is anyting other than 6 exit Matlab
if rank(ctrb(sys))~= 6
   exit(0);
end

% Weight matrix Q, R can be adjusted according to ones wish.
% Weights we put on the State X of the matrix
Q = 1*eye(6);

% Weights we put on the Control matrix [Ux,Uy, Uz]. Here Ux, Uy and Uz are
% the accelerations in the respective direction which will drive the chaser satellite. 
R = 10^4*eye(3);

% Finding the Optimal Gain matrix K
[K,P,E] = lqr(A,B,Q,R);

% New System with Control Feedback
Ac = A-B*K;

% U = Control Variable at Various times
U = zeros(3,numpoints);

% X - State Vector
X(:,1) = [rho0_RSW;rhodot0_RSW];

% dt = time step for simulation 
dt = Ts;

for j = 1:dt:numpoints
   U(:,j) = -K*X(:,j);
   Xdot = A*X(:,j)+B*U(:,j);
   X(:,j+1) = X(:,j)+Xdot*dt;
   
   % Stop Condition
   if rms(X(1:3,j)) <= 0.0005
       endtime = j;
       break;
   else
       endtime = j;
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
p1 = plot(t_span(1:endtime),X(1,1:endtime),'r','LineWidth',2);
legend('LF')
xlabel('t (in s)');
ylabel('x (in Km)');
hold off

subplot(3,1,2)
hold on
grid on
p1 = plot(t_span(1:endtime),X(2,1:endtime),'r','LineWidth',2);
xlabel('t (in s)');
ylabel('y (in Km)');
hold off

subplot(3,1,3)
hold on
grid on
p1 = plot(t_span(1:endtime),X(3,1:endtime),'r','LineWidth',2);
xlabel('t (in s)');
ylabel('z (in Km)');
hold off

figure;

subplot(3,1,1)
hold on
grid on
p1 = plot(t_span(1:endtime),U(1,1:endtime),'r','LineWidth',2);
legend('LF')
xlabel('t (in s)');
ylabel('Ux (in KN)');
hold off

subplot(3,1,2)
hold on
grid on
p1 = plot(t_span(1:endtime),U(2,1:endtime),'r','LineWidth',2);
xlabel('t (in s)');
ylabel('Uy (in KN)');
hold off

subplot(3,1,3)
hold on
grid on
p1 = plot(t_span(1:endtime),U(3,1:endtime),'r','LineWidth',2);
xlabel('t (in s)');
ylabel('Uz (in KN)');
hold off