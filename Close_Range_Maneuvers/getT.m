function T = getT(f,e,h)

global mu;
T = [(1+e*cos(f))*eye(3) 0*eye(3);...
    -e*sin(f)*eye(3) h/(mu^2*(1+e*cos(f)))*eye(3)];

end