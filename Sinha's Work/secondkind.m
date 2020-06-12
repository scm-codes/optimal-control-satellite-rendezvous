function soln = secondkind(problem)
% soln = secondkind(problem)


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

U_previous = 1;
U_current = 4*(t-ti)/beta-2;
U_array = [1,4*(t-ti)/beta-2];

for i = 3:1:m
    U_next = 2*(2*(t-ti)/beta-1)*U_current-U_previous;
    U_array = [U_array,U_next];
    U_previous = U_current;
    U_current = U_next;
end
    
State = Kronecker_Product(eye(n),U_array)*expansioncoeffs_State;
State = vpa(expand(State),8);
soln.StateExpression = State;
soln.StateFunc = matlabFunction(State);
soln.StateExpnsCoeffs = expansioncoeffs_State;
    
Phi = Kronecker_Product(eye(n),U_array)*expansioncoeffs_Phi;
Phi = expand(Phi);
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
    Q = OpMat_Product(m,n,periodic_mat,ti,tf,beta);

    Z = beta*Kronecker_Product(aperiodic_mat,G') + beta*Kronecker_Product(eye(n),G')*Q;
    expansion_coeffs_X = (eye(size(aperiodic_mat).*size(G)) - Z)\State0;
    expansion_coeffs_Phi = (eye(size(aperiodic_mat).*size(G)) - Z)\Phi0;
end


%%%%%%%%%%% OpMat_Integration.m

function G = OpMat_Integration(m)
    % This code calculates the operational matrix of Integration
    % m = No. of terms in Chebyshev Expansion
    
    % U = Shifted Chebyshev Polynomials of Second Kind
    U = sym('U',[m 1]);
    
    % System of Equations to get Operational Matrix of Integration of
    % Shifted Chebyshev Polynomials
    eqn = [0.5*U(1) + 0.25*U(2),-3/8*U(1) + 1/8*U(3)]; 
    for i = 3:1:m
        if i == m 
            nxteqn = - 0.25*U(i-1)/i +(-1)^(i-1)*U(1)/(2*i);
            eqn = [eqn,nxteqn]; %#ok<AGROW>
        else  
            nxteqn = 0.25*U(i+1)/i - 0.25*U(i-1)/i +(-1)^(i-1)*U(1)/(2*i);
            eqn = [eqn,nxteqn]; %#ok<AGROW>
        end
    end
   
    % G - Operation Matrix of Integration
    [G,~] = equationsToMatrix(eqn,U(1:m));
end


%%%%%%%%%%%% OpMat_Product.m

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