function [PhiT,time] = getPhiT_firstkind(m,n,aperiodic_mat,periodic_mat,ti,tf,G)
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

T_previous = 1;
T_current = 2*(t-ti)/beta-1;
T_array = [1,2*(t-ti)/beta-1];

for i = 3:1:m
    T_next = 2*(2*(t-ti)/beta-1)*T_current-T_previous;
    T_array = [T_array,T_next];
    T_previous = T_current;
    T_current = T_next;
end

%%%% Exapnsion Coefficients
expansioncoeffs_Phi = getExpansion_Coeffs(m,n,aperiodic_mat,periodic_mat,Phi0,ti,tf,beta,G);

Phi = kron(eye(n),T_array)*expansioncoeffs_Phi;
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
    Q = OpMat_Product(m,n,periodic_mat,ti,tf);

    Z = beta*kron(aperiodic_mat,G') + beta*kron(eye(n),G')*Q;
    expansion_coeffs_Phi = (eye(size(aperiodic_mat).*size(G)) - Z)\Phi0;
end

%%%%%%%%%%%% OpMat_product
function Q = OpMat_Product(m,n,periodic_func,t_initial,t_final)
    % This code Calculates the Product Operational Matrix 
    % m = No. of terms in Chebyshev Expansion
    
    % Chebyshev coefficeeint Matrix
    Cheb_Coeff_Mat = ShiftedChebCoeff(2*m,n,periodic_func,t_initial,t_final);
    
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
                            q(j,i) = Cheb_Coeff_Vec(i,1)/2;
                        elseif i == j
                            q(j,i) = Cheb_Coeff_Vec(1,1) + Cheb_Coeff_Vec(2*(j-1)+1,1)/2;
                        else 
                            q(j,i) = 0.5*(Cheb_Coeff_Vec(i+j-1,1) + Cheb_Coeff_Vec(abs(i-j)+1,1));
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
function Coeff_Matrix = ShiftedChebCoeff(m,n,periodic_func,ti,tf)
    % This code computes the Shifted Chebyshev expansion coeffs for the
    % periodic terms in the system Matrix
    
    beta = tf - ti;
    % Coeffcient Matrix
    Coeff_Matrix = zeros(n*m,n);
    
    % Generating Coeff Shifted Chebyshev Polynomial
    for l = 1:1:n
        for k = 1:1:n


            for i = 0:1:m-1
    
                % Chebyshev Polynomial
                if i == 0
                    T_poly = @(y) 1;
                elseif i == 1
                    T_poly = @(y) 2*y-1;
                else 
                    T_previous = @(y) 1;
                    T_current = @(y) 2*y-1;

                    for j = 2:1:i
                        T_poly = @(y) 2*(2*y-1).*T_current(y) - T_previous(y);
                        T_previous = @(y) T_current(y);
                        T_current = @(y) T_poly(y);
                    end
                end

                % Finding Coefficient
                weight = @(y) 1./sqrt(y-y.^2);
                int_fun = @(y) weight(y).*periodic_func{k,l}(beta*y+ti).*T_poly(y);

               if i == 0
                   Coeff_Matrix((k-1)*m+i+1,l) = 1/pi*integral(int_fun,0,1,'RelTol',1e-9,'AbsTol',1e-9);
               else 
                   Coeff_Matrix((k-1)*m+i+1,l) = 1/(pi/2)*integral(int_fun,0,1,'RelTol',1e-9,'AbsTol',1e-9);
               end
            end
        end
    end
end