function uSoln= Fourier_ControlInterpolant(problem,params,tgrid,Zsoln)
    
    size_t = size(tgrid);
    if size_t(1) == 1
        tgrid = tgrid';
    end
    % For Code readability
    Fou = problem.fourier;
    
    % Unpack problem
    m = Fou.harmonics;
    n = Fou.frequency;
        
    % Set up Harmonics array
    j = (1:1:m);
    
    % trgtdata construction
    mu = 3.986004415000000e+05;
    params.target.tspan = tgrid;
    e = params.target.e;
    h = params.target.h;
    [tet,dtet,ddtet] = target_info(params.target);
    tet = tet';
    dtet = dtet';
    ddtet = ddtet';
    mt = params.chaser.mass;
    rt = h^2./(mu.*(1+e.*cos(tet)));
    
    
    % trigdata construction
    c = cos(n.*kron(j,tgrid)');
    s = sin(n.*kron(j,tgrid)');
    dc = -j'.*n.*s;
    ds = j'.*n.*c;
    ddc = -j'.*n.*ds;
    dds = j'.*n.*dc;
    
    %%%% Reconstruct the solution
    % Fourier Coefficients
    gam1 = Zsoln(6*m+1,:);
    gam2 = Zsoln(6*m+2,:);
    gam3 = Zsoln(6*m+3,:);
    
    alp1 = Zsoln(1:m,:);
    bet1 = Zsoln(m+1:2*m,:);
    alp2 = Zsoln(2*m+1:3*m,:);
    bet2 = Zsoln(3*m+1:4*m,:);
    alp3 = Zsoln(4*m+1:5*m,:);
    bet3 = Zsoln(5*m+1:6*m,:);

    % path
    q1 = gam1 + sum(alp1.*c)' + sum(bet1.*s)';
    q2 = gam2 + sum(alp2.*c)' + sum(bet2.*s)';
    q3 = gam3 + sum(alp3.*c)' + sum(bet3.*s)';
    
    dq1 = sum(alp1.*dc)' + sum(bet1.*ds)';
    dq2 = sum(alp2.*dc)' + sum(bet2.*ds)';
    dq3 = sum(alp3.*dc)' + sum(bet3.*ds)';
    
    ddq1 = sum(alp1.*ddc)' + sum(bet1.*dds)';
    ddq2 = sum(alp2.*ddc)' + sum(bet2.*dds)';
    ddq3 = sum(alp3.*ddc)' + sum(bet3.*dds)';
    
    % Non Flat Frame
    rc = sqrt((rt+q1).^2 + q2.^2 + q3.^2);
    cnst = mu./rc.^3;
    ux = mt.*(ddq1 - 2.*dtet.*dq2 - ddtet.*q2 - dtet.^2.*q1 - mu./rt.^2 + cnst.*(rt+q1));
    uy = mt.*(ddq2 + 2.*dtet.*dq1 + ddtet.*q1 - dtet.^2.*q2 + cnst.*q2);
    uz = mt.*(ddq3 + cnst.*q3);
    uSoln = [ux,uy,uz];
end