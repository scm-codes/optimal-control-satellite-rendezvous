function soln = firstkind(problem)
% soln = firstkind(problem)


%%%% Unpack problem
m = problem.m;
n = problem.n;
aperiodic_mat = problem.aperiodic_mat;
periodic_mat = problem.periodic_mat;
initialcond = problem.initialcond;
ti = problem.ti;
tf = problem.tf;

% beta
beta = tf-ti;

%%%% Initial Condition Set-Up
% State
State0 = zeros(m*n,1);
for i = 1:1:n
    State0(m*(i-1)+1,1) = initialcond(i,1);
end
% Phi
Phi0 = zeros(n*m,n);
for i = 1:1:n
    Phi0(m*(i-1)+1,i) = 1;
 end

%%%% Exapnsion Coefficients
[expansioncoeffs_State,...
    expansioncoeffs_Phi] = getExpansion_Coeffs(m,n,aperiodic_mat,periodic_mat,State0,Phi0,ti,tf,beta);

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
    
State = Kronecker_Product(eye(n),T_array)*expansioncoeffs_State;
State = expand(State);
soln.StateExpression = State;
soln.StateFunc = matlabFunction(State);
soln.StateExpnsCoeffs = expansioncoeffs_State;
    
Phi = Kronecker_Product(eye(n),T_array)*expansioncoeffs_Phi;
Phi = vpa(expand(Phi),8);
soln.PhiExpression = Phi;
soln.PhiFunc = matlabFunction(Phi);
soln.PhiExpnCoeffs = expansioncoeffs_Phi;

end


%-------------------------------------------------------------------------%
%                           Utillity Functions                            %
%-------------------------------------------------------------------------%

%%%%%%%%%% getExpansion_Coeffs.m

function [expansion_coeffs_X,expansion_coeffs_Phi] = getExpansion_Coeffs(m,n,aperiodic_mat, periodic_mat,State0,Phi0,ti,tf,beta)
    % Return the Exapnsion coeffs for the solution of State and Phi
    G = double(OpMat_Integration(m)); 
    Q = OpMat_Product(m,n,periodic_mat,ti,tf);

    Z = beta*Kronecker_Product(aperiodic_mat,G') + beta*Kronecker_Product(eye(n),G')*Q;
    expansion_coeffs_X = (eye(size(aperiodic_mat).*size(G)) - Z)\State0;
    expansion_coeffs_Phi = (eye(size(aperiodic_mat).*size(G)) - Z)\Phi0;
end


%%%%%%%%%%% OpMat_Integration.m

function G = OpMat_Integration(m)
    % This code calculates the operational matrix of Integration
    % m = No. of terms in Chebyshev Expansion
    
    % T = Shifted Chebyshev Polynomials of First Kind
    T = sym('T',[m 1]);
    
    % System of Equations to get Operational Matrix of Integration of
    % Shifted Chebyshev Polynomials
    eqn = [0.5*T(1) + 0.5*T(2),-1/8*T(1) + 1/8*T(3)]; 
    for i = 3:1:m
        if i == m 
            nxteqn = - 0.25*T(i-1)/(i-2) -(-1)^(i-1)*T(1)/(2*((i-1)^2 - 1));
            eqn = [eqn,nxteqn]; %#ok<AGROW>
        else  
            nxteqn = 0.25*T(i+1)/i - 0.25*T(i-1)/(i-2) -(-1)^(i-1)*T(1)/(2*((i-1)^2 - 1));
            eqn = [eqn,nxteqn]; %#ok<AGROW>
        end
    end
   
    % G - Operation Matrix of Integration
    [G,~] = equationsToMatrix(eqn,T(1:m));
end


%%%%%%%%%%%% OpMat_Product.m

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

function Coeff_Matrix = ShiftedChebCoeff(m,n,periodic_func,t_initial,t_final)
    % This code computes the Shifted Chebyshev expansion coeffs for the
    % periodic terms in the system Matrix
    
    beta = t_final - t_initial;
    
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
                int_fun = @(y) weight(y).*periodic_func{k,l}(beta*y+t_initial).*T_poly(y);

               if i == 0
                   Coeff_Matrix((k-1)*m+i+1,l) = 1/pi*integral(int_fun,0,1,'RelTol',1e-9,'AbsTol',1e-9);
               else 
                   Coeff_Matrix((k-1)*m+i+1,l) = 1/(pi/2)*integral(int_fun,0,1,'RelTol',1e-9,'AbsTol',1e-9);
               end
            end
        end
    end
end


%%%%%%%%%%%%%% Kronecker_Product.m

function Kronecker_Mat = Kronecker_Product(A,B)
    % Returns the Kronecker product of Square Matrices A and B
    n = length(A);
    
    Kronecker_Mat = [];
    for i = 1:1:n
        Mat = [];
        for j = 1:1:n
            Mat = [Mat, A(i,j)*B]; %#ok<AGROW>
        end
        Kronecker_Mat = [Kronecker_Mat;Mat];  %#ok<AGROW>
    end
end