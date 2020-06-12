function [c] = InitialClassicalOrbitalElements(satellite)
% Sai Charan Malladi
% AE16B029 IIT Madras
% 2/1/2019
% This function sets the initial conditions in terms of classical orbital
% elements.

% Please enter the following in the given units only
% a = semi major axis of orbit [km]
% e = eccentricity 
% i = inclination of the orbit [rad]
% raan = Right ascesion of ascending node [rad]
% ap = argument of perigee [rad]
% theta = true anomaly [rad]

switch satellite
    
        %%%% Case 1    
    case 'ChaserCase1'
        c.a = 6835;
        c.e = 0;
        c.i = 0;
        c.raan = 0;
        c.ap = 0;
        c.theta = 0;
    case'TargetCase1'
        c.a = 6815;
        c.e = 0;
        c.i = 0;
        c.raan = 0;
        c.ap = 0;
        c.theta = 0;
    
        %%%% Case 2
    case 'ChaserCase2'
        c.a = 11000;
        c.e = 0.111;
        c.i = 0.01*pi/180;
        c.raan = 0;
        c.ap = 0;
        c.theta = 0;
    case 'TargetCase2' 
        c.a = 11010;
        c.e = 0.1125;
        c.i = 0;
        c.raan = 0;
        c.ap = 0;
        c.theta = 0;
        
    case 'ChaserCase3'
        c.a = 7440;
        c.e = 0.1;
        c.i = 0;
        c.raan = 0;
        c.ap = 0;
        c.theta = 0;
    case 'TargetCase3' 
        c.a = 7420;
        c.e = 0.12;
        c.i = 0;
        c.raan = 0;
        c.ap = 0;
        c.theta = 0;
end 