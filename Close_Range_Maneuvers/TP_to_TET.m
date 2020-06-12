function [theta,M] = TP_to_TET(t,e,a)
% Sai Charan Malladi
% AE16B029 IIT Madras
% 2/1/2019
% This function converts time since periapsis(t) to true anomaly(theta)

global mu;

% Calculation of Mean motion n
n = sqrt(mu/a^3);

%Calulation of Mean anomaly M
M = n*(t);

% Mean Anomaly to Eccentric Anomaly
M2E  = @(M,e) fzero(@(E) E-e*sin(E)-M, M);

E = wrapTo2Pi(M2E(M,e));

%Calculatiom of theta
theta = wrapTo2Pi(2*atan(sqrt((1+e)/(1-e))*tan(E/2)));
end