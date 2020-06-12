function [thetadot,thetadotdot] = DERIVATIVES_of_THETA(theta,e,a,raan,i,ap,h)
% Sai Charan Malladi
% AE16B029 IIT Madras
% 2/1/2019
% This function computes the derivatives of true anomaly

% Compute r and v vectors in ECI frame at a given theta 
[r,v] = SV_from_COE(h,e,raan,i,ap,theta);

% Compute r.r
modr = norm(r);

% Derivatives
thetadot = h/modr^2;
thetadotdot = -2*dot(v,r)*thetadot/modr^2;
end