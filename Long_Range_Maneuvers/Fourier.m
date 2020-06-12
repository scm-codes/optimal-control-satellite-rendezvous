function soln = Fourier(problem,params)

    % For Code readability
    Opt = problem.options;
    Fou = problem.fourier;
    B = problem.bounds;
    T = problem.time;
    
    % Unpack problem
    tSpan = T.limits;
    tSegments = T.segments;
    m = Fou.harmonics;
    n = Fou.frequency;
    Jval = Fou.Objective;
    tLen = tSegments+1;
    
    % Set up time for optimzation
     tgrid = (linspace(tSpan(1),tSpan(2),tLen))';
        
    % Set up Harmonics array
    j = (1:1:m);
    
    % trgtdata construction
    global mu;
    params.target.tspan = tgrid;
    e = params.target.e;
    h = params.target.h;
    [tet,dtet,ddtet] = target_info(params.target);
    tet = tet';
    dtet = dtet';
    ddtet = ddtet';
    mt = params.chaser.mass;
    rt = h^2./(mu.*(1+e.*cos(tet)));
    
    % Unpack trgtdata
    trgtdata.mt = mt;
    trgtdata.rt = rt;
    trgtdata.dtet = dtet;
    trgtdata.ddtet = ddtet;
    
    % trigdata construction
    c = cos(n.*kron(j,tgrid)');
    s = sin(n.*kron(j,tgrid)');
    dc = -j'.*n.*s;
    ds = j'.*n.*c;
    ddc = -j'.*n.*ds;
    dds = j'.*n.*dc;
    
    trigdata.c = c;
    trigdata.s = s;
    trigdata.dc = dc;
    trigdata.ds = ds;
    trigdata.ddc = ddc;
    trigdata.dds = dds;
    
    % Guess trjectory
    if ~isfield(problem, 'guess')
        Zguess = [zeros(6*m,1);166.4;0;0];
    else
        G = problem.guess;
        Zguess = interp1(G.time,G.fouriercoeffs,tgrid);
    end
    
    %%%% Set up the problem for fmincon
    
    % Objective Functio
    P.objective = @(Z)(...
            PathObj_Fourier_Flat(Z,m,j,Jval));
    
    % Constraits Function
    P.nonlcon = @(Z)(...
        Cst_Fourier(Z,m,tLen,trigdata,trgtdata,B));
    
    P.x0 = Zguess;
    P.nvars = 6*m+3;
    P.x0 = Zguess;
    P.lb = [];
    P.ub = [];
    P.Aineq = []; P.bineq = [];
    P.Aeq = []; P.beq = [];
    P.options = Opt.nlpOpt;
    P.solver = 'fmincon';

    %%%% Call fmincon to solve the non-linear program (NLP)
    tic;
    [Zsoln, ~,exitFlag,output] = fmincon(P);
    nlpTime = toc;
    
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
    
    xSoln = [q1,q2,q3,dq1,dq2,dq3,ddq1,ddq2,ddq3];
    
    % Non Flat Frame
    rc = sqrt((rt+q1).^2 + q2.^2 + q3.^2);
    cnst = mu./rc.^3;
    ux = mt.*(ddq1 - 2.*dtet.*dq2 - ddtet.*q2 - dtet.^2.*q1 - mu./rt.^2 + cnst.*(rt+q1));
    uy = mt.*(ddq2 + 2.*dtet.*dq1 + ddtet.*q1 - dtet.^2.*q2 + cnst.*q2);
    uz = mt.*(ddq3 + cnst.*q3);
    uSoln = [ux,uy,uz];
    U = (ux.^2 + uy.^2 + uz.^2);
    Obj = trapz(tgrid,U);
   
    %%%% Construct the interpolation
    %%%% Store the results
    soln.grid.time = tgrid;
    soln.grid.fouriercoeffs = Zsoln;
    soln.grid.state = xSoln;
    soln.grid.control = uSoln;
    soln.info = output;
    soln.info.nlpTime = nlpTime;
    soln.info.exitFlag = exitFlag;
    soln.info.objVal = Obj;
    soln.problem = problem; % Return the fully detailed problem struct
    
end