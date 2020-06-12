function [rho_RSW, rhodot_RSW] = RSW_from_SV(r_ECI, v_ECI,rho_ECI,rhodot_ECI)
% Sai Charan Malladi
% AE16B029 IIT Madras
% 2/1/2019
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