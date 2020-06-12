function dX = LERM(t,X,tp,e,a,raan,i,ap,h)
% Sai Charan Malladi
% AE16B029 IIT Madras
% 2/1/2019
% This function is for ODE45 solver to solve the linear equations of relative motion
global mu;

%Compute theta and its derivative at each time as the coefficients of the
%linear equations vary with time.
[theta,~] = TP_to_TET(t.*1+tp,e,a);
r = h^2/mu/(1+e*cos(theta));
[thetadot,thetadotdot] = DERIVATIVES_of_THETA(theta,e,a,raan,i,ap,h);

%Position
dX(1,1) = X(4);
dX(2,1) = X(5);
dX(3,1) = X(6);

%Velocity
dX(4,1) = (thetadot^2 + 2*mu/r^3)*X(1) + thetadotdot*X(2) + 2*thetadot*X(5);
dX(5,1) = -thetadotdot*X(1) + (thetadot^2 - mu/r^3)*X(2) - 2*thetadot*X(4);
dX(6,1) = -mu/r^3*X(3);
end