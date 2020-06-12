function [Phi,R_bar,U_min,U_max,C_m,Rs,F,Z1,Z2] = getMPCstuff(N_p,N_c,A_cm,B_cm,C_cm,D_cm,T,Xsetpoint,W,Umin,Umax)
    % Fuction prepares all relevant matrices for optimization
    
    sizeAcm = size(A_cm);
    % Discrete state-space representation
    A_m = expm(A_cm*T);
    B_m = integral(@(x) expm(A_cm*x)*B_cm,0,T,'ArrayValued',true);
    C_m = C_cm;
    D_m = D_cm;
    
    % Increment Input output model(IIO)
    A = [   A_m        zeros(sizeAcm);
            C_m*A_m    eye(sizeAcm)  ];
    B = [   B_m;
            C_m*B_m     ];
    C = [zeros(sizeAcm) eye(sizeAcm)];
         
    % Calculation of F and Phi matrices
    F = C*A;
    for i = 2:N_p
       F = [F;C*A^i]; 
    end
    
    sizeCB = size(C*B);
    Phi = [];
    for i = 1:N_p
        phi = [];
        for j = 1:N_c
            if j > i
                phi = [phi zeros(sizeCB)];
            else
                phi = [phi C*A^(i-j)*B];
            end
        end
        Phi = [Phi;phi];
    end
    
    % Vector that stores set point information Rs
    Rs_bar = [];
    for i = 1:N_p
        Rs_bar = [Rs_bar eye(sizeAcm)];
    end
    Rs = Rs_bar'*Xsetpoint;
    
    % Vector that puts weight on control R_bar
    sizeW = size(W);
    R_bar = [];
    for i = 1:N_c
        r_bar = [];
        for j = 1:N_c
            if j == i
                r_bar = [r_bar W];
            else
                r_bar = [r_bar zeros(sizeW)];
            end
        end
        R_bar = [R_bar; r_bar];
    end
    
    %%%% Constraint Formulation
    Z1 = [];
    for i = 1:N_c
        Z1 = [Z1;eye(sizeW)];
    end
    
    Z2 = [];
    for i = 1:N_c
        z2 = [];
        for j = 1:N_c
            if j <= i
                z2 = [z2 eye(sizeW)];
            else z2 = [z2 zeros(sizeW)];
            end
        end
        Z2 = [Z2;z2];
    end
    
    % Concatenated coloumn block vectors
    U_min = [];
    U_max = [];
    
    for i = 1:N_c
        U_min = [U_min; Umin];
        U_max = [U_max; Umax];
    end
end