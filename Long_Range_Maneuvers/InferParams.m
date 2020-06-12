function [h,T,n,tp] = InferParams(p)

global mu;

% Unpack p
a = p.a;
e = p.e;
tet0 = p.theta0;

%%%% Angular momentum per unitmass
h = sqrt(a*mu*(1-e^2));

%%%% Time period
T = 2*pi*sqrt(a^3/mu);

%%%% Mean motion
n = 2*pi/T;

%%%% Time since periapsis
%Calculation of Eccentric Anomaly E
E0 = wrapTo2Pi(2*atan(sqrt((1-e)/(1+e)))*tan(tet0/2));

%Calculation of Mean Anomaly M
M0 = wrapTo2Pi(E0-e*sin(E0));

%Calculation of Mean Motion n
n = sqrt(mu/a^3);

tp = M0/n;


end