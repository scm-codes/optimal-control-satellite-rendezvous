function [theta,dtheta,ddtheta] = target_info(p)
% [theta,dtheta,ddtheta] = ThetaInfo(p)
% Utility function of RendezvousDynamics.m
%

global mu;

% unpack params
n = p.n;
e = p.e;
h = p.h;
i = p.i;
raan = p.raan;
ap = p.ap;
tp = p.tp;
t = p.tspan;

tnet = t+tp;
nt = length(t);

% Calulation of Mean anomaly M
M = n.*(tnet);

% Mean Anomaly to Eccentric Anomaly
M2E  = @(M,e) fzero(@(E) E-e.*sin(E)-M, M);
E = zeros(1,nt);
for j = 1:nt
    E(j) = wrapTo2Pi(M2E(M(j),e));
end

%%%% Calculatiom of theta
theta = wrapTo2Pi(2.*atan(sqrt((1+e)/(1-e)).*tan(E./2)));

% rp and vp are column vectors in perifocal frame
rp = (h^2/mu) * (1./(1 + e.*cos(theta))) .* (cos(theta).*[1;0;0] + sin(theta).*[0;1;0]); 
vp = (mu/h) * (-sin(theta).*[1;0;0] + (e + cos(theta)).*[0;1;0]);

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
r = zeros(3,nt);
v = zeros(3,nt);
for j = 1:nt
    r(:,j) = transmzat*rp(:,j);
    v(:,j) = transmzat*vp(:,j);
end

% Compute r.r
modr = zeros(1,nt);
for j = 1:nt
    modr(j) = norm(r(:,j));
end

%%%% Derivatives
dtheta = h./modr.^2;
ddtheta = -2.*dot(v,r).*dtheta./modr.^2;

end