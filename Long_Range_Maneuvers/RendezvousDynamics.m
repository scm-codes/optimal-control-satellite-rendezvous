function dX = RendezvousDynamics(t,X,U,params)
% dX = RendezvousDynamics(t,X,params)
%
%

global mu;

% unpack the state
x = X(1,:);
y = X(2,:);
z = X(3,:);
dx = X(4,:);
dy = X(5,:);
dz = X(6,:);

% unpack control
ux = U(1,:);
uy = U(2,:);
uz = U(3,:);

% unpack params
h = params.target.h;
e = params.target.e;
theta = params.target.theta;
dtheta = params.target.dtheta;
ddtheta = params.target.ddtheta;
tspan = params.target.tspan;
mc = params.chaser.mass;

% calculate theta and derivatives of target
tet = interp1(tspan,theta,t);
dtet = interp1(tspan,dtheta,t);
ddtet = interp1(tspan,ddtheta,t);

% distance of target from Earth cetered frame
r_t = h^2./(mu.*(1+e.*cos(tet)));

% distance of chaser from Earth centered frame
r_c = sqrt((r_t+x).^2 + y.^2 + z.^2);

cnst = 1./r_c.^3;
% Dynamics
ddx = (dtet.^2-mu.*cnst).*x + ddtet.*y + 2.*dtet.*dy....
    + mu.*(1./r_t.^2 - r_t.*cnst) + ux./mc;
ddy = ddtet.*x + (dtet.^2-mu.*cnst).*y + 2.*dtet.*dx + uy./mc;
ddz = -mu.*cnst.*z + uz./mc;

% Pack
dX = [dx;dy;dz;ddx;ddy;ddz];

end