function [c, ceq] = Cst_Fourier(Z,m,tLen,trigdata,trgtdata,B)
    % Constraint Function for Fourier Collocation
    
    % Unpack trgtdata
    mt = trgtdata.mt;
    rt = trgtdata.rt;
    dtet = trgtdata.dtet;
    ddtet = trgtdata.ddtet;
    
    % Unpack trigdata
    c = trigdata.c;
    s = trigdata.s;
    dc = trigdata.dc;
    ds = trigdata.ds;
    ddc = trigdata.ddc;
    dds = trigdata.dds;
    
    % Unpack B
    Xi = B.initialstate;
    Xf_low = B.finalstate.low;
    Xf_upp = B.finalstate.upp;
    X_low = B.state.low;
    X_upp = B.state.upp;
    U_low = B.control.low;
    U_upp = B.control.upp;
    
    % Fourier Coefficients
    gam1 = Z(6*m+1,:);
    gam2 = Z(6*m+2,:);
    gam3 = Z(6*m+3,:);
    
    alp1 = Z(1:m,:);
    bet1 = Z(m+1:2*m,:);
    alp2 = Z(2*m+1:3*m,:);
    bet2 = Z(3*m+1:4*m,:);
    alp3 = Z(4*m+1:5*m,:);
    bet3 = Z(5*m+1:6*m,:);

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
    global mu;
    rc = sqrt((rt+q1).^2 + q2.^2 + q3.^2);
    cnst = mu./rc.^3;
    ux = mt.*(ddq1 - 2.*dtet.*dq2 - ddtet.*q2 - dtet.^2.*q1 - mu./rt.^2 + cnst.*(rt+q1));
    uy = mt.*(ddq2 + 2.*dtet.*dq1 + ddtet.*q1 - dtet.^2.*q2 + cnst.*q2);
    uz = mt.*(ddq3 + cnst.*q3);
    
    
    %%%% Boundary Constraints Formulation
    % Initial
    cstXi_upp = [q1(1);q2(1);q3(1);dq1(1);dq2(1);dq3(1)] - Xi;
    cstXi_low = Xi - [q1(1);q2(1);q3(1);dq1(1);dq2(1);dq3(1)];
    
    % Final
    cstXf_upp = [q1(end);q2(end);q3(end);dq1(end);dq2(end);dq3(end)] - Xf_upp;
    cstXf_low = Xf_low - [q1(end);q2(end);q3(end);dq1(end);dq2(end);dq3(end)];
    c_bnd = [cstXi_upp;cstXi_low;cstXf_upp;cstXf_low];
    ceq_bnd = [];
    
    %%%% Path Constraints
    cstX_upp = [q1(2:end-1);q2(2:end-1);q3(2:end-1);dq1(2:end-1);dq2(2:end-1);dq3(2:end-1)] - kron(X_upp,ones(tLen-2,1));
    cstX_low = kron(X_low,ones(tLen-2,1)) - [q1(2:end-1);q2(2:end-1);q3(2:end-1);dq1(2:end-1);dq2(2:end-1);dq3(2:end-1)];
    cstU_upp = [ux(2:end-1);uy(2:end-1);uz(2:end-1)] - kron(U_upp.*ones(3,1),ones(tLen-2,1));
    cstU_low = kron(U_low.*ones(3,1),ones(tLen-2,1)) - [ux(2:end-1);uy(2:end-1);uz(2:end-1)];
    c_path = [cstX_upp;cstX_low;cstU_upp;cstU_low];
    ceq_path = [];
    %%%% Collect all constraints
    c = [c_bnd;c_path];
    ceq = [ceq_bnd;ceq_path];
end