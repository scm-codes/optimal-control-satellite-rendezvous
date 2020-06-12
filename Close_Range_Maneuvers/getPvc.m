function Pvc = getPvc(theta,M,thetadot,n)
% Sai Charan Malladi
% AE16B029 IIT Madras
% 2/1/2019
% Returns the Pvc matrix as per Eq. 3.4 in Sherrill's paper

% Calcuation of R_vc and Omega_vc/c matrices as per Eq. 3.3 in Sherrill's
% paper.

R_vc = [cos(theta - M) sin(theta - M) 0;...
      -sin(theta - M) cos(theta - M) 0;...
            0              0         1];
        
Omega_vc = [     0    (thetadot - n)    0;...
           (n - thetadot)    0          0;...
                0            0          0];
            
% Inverse of Pvc matrix as per Eq. 3.4 in Sherrilli's paper

InvPvc = [      transpose(R_vc),               zeros(3);...
          -Omega_vc*transpose(R_vc),      transpose(R_vc)];
      
% Pvc calculation
Pvc = InvPvc';
end