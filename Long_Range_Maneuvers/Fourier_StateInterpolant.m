function xSoln= Fourier_StateInterpolant(problem,tgrid,Zsoln)
    
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
    
    % trigdata construction
    c = cos(n.*kron(j,tgrid)');
    s = sin(n.*kron(j,tgrid)');
    dc = -j'.*n.*s;
    ds = j'.*n.*c;
    
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
    
    xSoln = [q1,q2,q3,dq1,dq2,dq3];
end