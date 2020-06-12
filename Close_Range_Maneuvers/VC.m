function [rho, rhodot] = VC(theta_range,M_range,rho0,rhodot0,n,thetadot_range,t_range)
% Sai Charan Malladi
% AE16B029 IIT Madras
% 2/1/2019
% This function solves the Virtual Cheif as outiled in Sherill, "Dynamics
% and Control of Satellite Relative Motion in Elliptic Orbits using
% Lyapunov-Floquet Theory." The equation numbers are reference to this document.

% Initial state vector and other variables
X0 = [rho0;rhodot0];
theta0 = theta_range(1);
M0 = M_range(1);
thetadot0 = thetadot_range(1);

% Get the initial rotation matrix
Pvc0 = getPvc(theta0,M0,thetadot0,n);

% Compute the state vector for every iteration

for i = 1:1:length(theta_range)
    
    theta = theta_range(i);
    M = M_range(i);
    thetadot = thetadot_range(i);
    t = t_range(i);
    
    % Get the rotational matrix at every time instant
    Pvc = getPvc(theta,M,thetadot,n);
    
    % Get the Matrix Exponential matexp
    matexp =  [4-3*cos(n*t)      0   0           1/n*sin(n*t)     2/n*(1-cos(n*t))       0
               6*(sin(n*t)-n*t)  1   0           2/n*(cos(n*t)-1) 1/n*(4*sin(n*t)-3*n*t) 0
               0                 0   cos(n*t)    0                0                      1/n*sin(n*t)
               3*n*sin(n*t)      0   0           cos(n*t)         2*sin(n*t)             0
               6*n*(cos(n*t)-1)  0   0           -2*sin(n*t)      (4*cos(n*t)-3)         0
               0                 0  -n*sin(n*t)  0                0                      cos(n*t)];
           
    % Find the state at new time as per Eq. 3.5 in Sherrill's Paper
    X = Pvc*matexp*Pvc0'*X0;
    
    % Split the state vector into rho and rhodot
    rho(:,i) = X(1:3);
    rhodot(:,i) = X(4:6);
end
end 