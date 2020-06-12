function tp = TET_to_TP(theta,e,a)
% Sai Charan Malladi
% AE16B029 IIT Madras
% 2/1/2019
% This function converts True anomaly(theta) to Time since periapsis(tp)
global mu;

%Calculation of Eccentric Anomaly E
E = wrapTo2Pi(2*atan(sqrt((1-e)/(1+e)))*tan(theta/2));

%Calculation of Mean Anomaly M
M = wrapTo2Pi(E-e*sin(E));

%Calculation of Mean Motion n
n = sqrt(mu/a^3);

tp = M/n;
end

