function [PhiT,time] = getPhiT_secondkind(m,n,aperiodic_mat,periodic_mat,ti,tf,G)
% [PhiT,time] = getPhiT_firstkind(m,n,aperiodic_mat,periodic_mat,ti,tf,G)

% beta
beta = tf-ti;

%%%% Initial Condition Set-Up
tic;

% Phi
Phi0 = zeros(n*m,n);
for i = 1:1:n
    Phi0(m*(i-1)+1,i) = 1;
end
 
%%%% 
syms t

U_previous = 1;
U_current = 4*(t-ti)/beta-2;
U_array = [1,4*(t-ti)/beta-2];

for i = 3:1:m
    U_next = 2*(2*(t-ti)/beta-1)*U_current-U_previous;
    U_array = [U_array,U_next];
    U_previous = U_current;
    U_current = U_next;
end

%%%% Exapnsion Coefficients
expansioncoeffs_Phi = getExpansion_Coeffs(m,n,aperiodic_mat,periodic_mat,Phi0,ti,tf,beta,G);

Phi = kron(eye(n),U_array)*expansioncoeffs_Phi;
PhiFunc = matlabFunction(Phi);
PhiT = PhiFunc(tf);
time = toc;
end

%-------------------------------------------------------------------------%
%                           Utillity Functions                            %
%-------------------------------------------------------------------------%

%%%%%%%%%% getExpansion_Coeffs.m
function expansion_coeffs_Phi = getExpansion_Coeffs(m,n,aperiodic_mat, periodic_mat,Phi0,ti,tf,beta,G)
    % Return the Exapnsion coeffs for the solution of State and Phi
    Q = OpMat_Product(m,n,periodic_mat,ti,tf,beta);

    Z = beta*kron(aperiodic_mat,G') + beta*kron(eye(n),G')*Q;
    expansion_coeffs_Phi = (eye(size(aperiodic_mat).*size(G)) - Z)\Phi0;
end

%%%%%%%%%%%% OpMat_product
function Q = OpMat_Product(m,n,periodic_func,ti,tf,beta)
    % This code Calculates the Product Operational Matrix 
    % m = No. of terms in Chebyshev Expansion
    
    % Chebyshev coefficeeint Matrix
    Cheb_Coeff_Mat = ShiftedChebCoeff(2*m,n,periodic_func,ti,tf,beta);
    
    % Q - Product Operation Matrix
    Q = zeros(n*m,n*m);
    
    for k = 1:1:n
        for l = 1:1:n
            Cheb_Coeff_Vec = Cheb_Coeff_Mat(2*m*(l-1)+1:2*m*l,k);
            q = zeros(m,m);
            
            % sub-matrices construction
            for i = 1:1:m
                for j = 1:1:m
                    if i == 1
                        q(j,i) = Cheb_Coeff_Vec(j,1);
                    else
                        if j == 1
                            q(j,i) = Cheb_Coeff_Vec(i,1);
                        elseif i == j
                            q(j,i) = q(abs(j-2)+1,abs(j-2)+1) + Cheb_Coeff_Vec(2*(j-1)+1,1);
                        else 
                            q(j,i) = q(abs(j-2)+1,abs(i-2)+1) + Cheb_Coeff_Vec(abs(i+j)-1,1);
                        end
                    end
                end
            end
            
            %Stacking sub-matrices
            Q(m*(l-1)+1:m*l,m*(k-1)+1:m*k) = q;
            
        end  
    end
end

%%%%%%%%%%%% ShiftedChebCoeff.m
function Coeff_Matrix = ShiftedChebCoeff(m,n,periodic_func,ti,tf,beta)
    % This code computes the Shifted Chebyshev expansion coeffs for the
    % periodic terms in the system Matrix
    
    weight = @(y) sqrt(y-y.^2);

U_cell = cell(2*m,1);

U_cell{1,1} = @(y) 1;
U_cell{2,1} = @(y) 4*y-2;

for i = 3:1:2*m
    U_cell{i,1} = @(y) 2*(2*y-1).*U_cell{i-1,1}(y)-U_cell{i-2,1}(y);
end

    % Coeffcient Matrix
    Coeff_Matrix = zeros(n*m,n);
    
    % Generating Coeff Shifted Chebyshev Polynomial
    for l = 1:1:n
        for k = 1:1:n
            for i = 0:1:m-1
    
                % Chebyshev Polynomial of Second Kind
                if i == 0
                    U_poly = @(y) 1;
                elseif i == 1
                    U_poly = @(y) 4*y-2;
                else
                    U_previous = @(y) 1;
                    U_current = @(y) 4*y-2;

                    for j = 2:1:i
                        U_poly = @(y) 2*(2*y-1).*U_current(y) - U_previous(y);
                        U_previous = @(y) U_current(y);
                        U_current = @(y) U_poly(y);
                    end
                end

                % Finding Coefficient
                weight = @(y) sqrt(y-y.^2);
                int_fun = @(y) weight(y).*periodic_func{k,l}(beta*y+ti).*U_poly(y);

                Coeff_Matrix((k-1)*m+i+1,l) = 8/pi*integral(int_fun,0,1,'RelTol',1e-9,'AbsTol',1e-9);
            end
        end
    end
end