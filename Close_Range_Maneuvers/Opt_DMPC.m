
function U = Opt_DMPC(Xk,Xk1,Uk1,Phi,R_bar,U_min,U_max,W,C_m,Rs,F,Z1,Z2)
    % Optimization function for Discrete MPC 
    
    % Augmented State matrix
    Xk_aug = [Xk-Xk1;C_m*Xk];
    
    sizeW = size(W);
    %%%% Cost function formulation
    H = Phi'*Phi + R_bar;
    f = -Phi'*(Rs - F*Xk_aug);

    % Mconst and Nconst matrices Mconst*delU <= Nconst 
    Mconst = [-Z2; Z2];
    Nconst = [-U_min + Z1*Uk1; U_max - Z1*Uk1];
    
    %%%% Constrained Optimization
    [delU,~] = quadprog(H,f,Mconst,Nconst);
    
    U = delU(1:sizeW(1),1);
    
end

