function X0 = rho0_info(p)
% X0 = rho0_info(p)
% 

% Initial State of target and chaser
[r0_t,v0_t] = SV_from_COE(p.target);
[r0_c,v0_c] = SV_from_COE(p.chaser);

% relative displacement vector in Earth-centered frame
rho0_ECI = r0_c - r0_t;
drho0_ECI = v0_c - v0_t;

% relative displacement vector in target frame
[rho0_RSW,drho0_RSW] = RSW_from_SV(r0_t,v0_t,rho0_ECI,drho0_ECI);
X0 = [rho0_RSW;drho0_RSW];

end

function [r, v] = SV_from_COE(p)
% [r, v] = SV_from_COE(p)
% This function computes the state vector (r,v) from the
% classical orbital elements (coe). 

global mu;

% unpack params
e = p.e;
h = p.h;
i = p.i;
raan = p.raan;
ap = p.ap;
tet0 = p.theta0;

% rp and vp are column vectors in perifocal frame
rp = (h^2/mu) * (1/(1 + e*cos(tet0))) * (cos(tet0)*[1;0;0] + sin(tet0)*[0;1;0]);
vp = (mu/h) * (-sin(tet0)*[1;0;0] + (e + cos(tet0))*[0;1;0]);

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

function [rho_RSW, rhodot_RSW] = RSW_from_SV(r_ECI, v_ECI,rho_ECI,rhodot_ECI)
% [rho_RSW, rhodot_RSW] = RSW_from_SV(r_ECI, v_ECI,rho_ECI,rhodot_ECI)
% This function converts states from ECI to RSW frame

%each of the components must be unit vectors
%radial componet
rvec = r_ECI./norm(r_ECI);

% cross-track component
wvec = cross(r_ECI, v_ECI);
wvec = wvec./norm(wvec);

% along-track component
svec = cross(wvec,rvec);
svec = svec./norm(svec);

% Transformation Matrix
transmat = [rvec'; svec'; wvec'];

% Calculate the instantaneous angular velocity
angvel_RSW = [0;0;norm(cross(r_ECI,v_ECI))/(norm(r_ECI))^2];

% relative position in RSW frame
rho_RSW = transmat*rho_ECI;
rhodot_RSW = transmat*rhodot_ECI - cross(angvel_RSW,rho_RSW);

end