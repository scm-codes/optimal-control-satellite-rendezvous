function [r, v] = SV_from_COE(h,e,raan,i,ap,theta)
% Sai Charan Malladi
% AE16B029 IIT Madras
% 2/1/2019
% This function computes the state vector (r,v) from the
% classical orbital elements (coe). 

global mu;

% rp and vp are column vectors in perifocal frame
rp = (h^2/mu) * (1/(1 + e*cos(theta))) * (cos(theta)*[1;0;0] + sin(theta)*[0;1;0]);
vp = (mu/h) * (-sin(theta)*[1;0;0] + (e + cos(theta))*[0;1;0]);

% Rotation matrix about the z-axis through the angle w
R3_W = [ cos(raan) sin(raan) 0
        -sin(raan) cos(raan) 0
        0 0 1];
% Rotation matrix about the x-axis through the angle i
R1_i = [1 0 0
        0 cos(i) sin(i)
        0 -sin(i) cos(i)];
    
% Rotation about the z-axis through the angle ap
R3_w = [ cos(ap) sin(ap) 0
        -sin(ap) cos(ap) 0
         0 0 1];
% Transformation matrix from COE to ECI
transmzat = R3_W'*R1_i'*R3_w';

% Compute r, v in ECI
r = transmzat*rp;
v = transmzat*vp;
end