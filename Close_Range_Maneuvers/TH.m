function [rho, rhodot] = TH(rho0,rhodot0,theta_range,a,e,M_range,h,n)
% Sai Charan Malladi
% AE16B029 IIT Madras
% 2/1/2019
% This function solves the Tschauner-Hempel Equations as outiled in
% Sherill, "Dynamics and Control of Satellite Relative Motion
% in Elliptic Orbits using Lyapunov-Floquet Theory." The equation numbers
% are reference to this document.

% format state vector and other variables
X0 = [rho0;rhodot0];
theta0 = theta_range(1);
M0 = M_range(1);

% find orbit parameter
p = a*(1-e^2);

% find constants and rotation matrices per Eq. 2.42 and Eq. 2.46 and 2.49
K0 = M0/(1-e^2)^(3/2);
PSI0 = getPSI(K0,theta0,e);
T0 = [(1+e*cos(theta0))*eye(3),             zeros(3);
    -e*sin(theta0)*eye(3),      p^2/(h*(1+e*cos(theta0)))*eye(3)];

% initialize variables
rho = nan(3,length(theta_range));
rhodot = nan(3,length(theta_range));

for i = 1:length(theta_range)
theta = theta_range(i);
M = M_range(i);

% get K per Eq. 2.49
K = M/(1-e^2)^(3/2);
    
% get PSI per Eq. 2.46
PSI = getPSI(K,theta,e);

% find T from Eq. 2.43
T = [(1+e*cos(theta))*eye(3),            zeros(3);
    -e*sin(theta)*eye(3),       p^2/(h*(1+e*cos(theta)))*eye(3)];

% find state at new time Eq. 2.51
X = inv(T)*PSI*inv(PSI0)*T0*X0;

% split into rho and rhodot
rho(:,i) = X(1:3);
rhodot(:,i) = X(4:6);
end


end

