function [rf,vf] = HCW(r0,v0,n,t)
% Sai Charan Malladi
% AE16B029 IIT Madras
% 2/1/2019
% Save Components to Variables
x0    = r0(1);
y0    = r0(2);
z0    = r0(3);
xdot0 = v0(1);
ydot0 = v0(2);
zdot0 = v0(3);

% The below equations are the explicit form of solutions for HCW model.

% Position
xf = 4*x0+2/n*ydot0+xdot0/n*sin(n*t)-(3*x0+2/n*ydot0)*cos(n*t);
yf = (6*x0+4*ydot0/n)*sin(n*t)+2*xdot0/n*cos(n*t)-(6*n*x0+3*ydot0)*t+...
     (y0-2*xdot0/n);
zf = z0*cos(n*t)+zdot0/n*sin(n*t);

% Velocity
xdotf = xdot0*cos(n*t)+(3*n*x0+2*ydot0)*sin(n*t);
ydotf = (6*n*x0+4*ydot0)*cos(n*t)-2*xdot0*sin(n*t)-(6*n*x0+3*ydot0);
zdotf = -z0*n*sin(n*t)+zdot0*cos(n*t);

% Save to Output Vectors
rf = [xf;yf;zf];
vf = [xdotf;ydotf;zdotf];
end