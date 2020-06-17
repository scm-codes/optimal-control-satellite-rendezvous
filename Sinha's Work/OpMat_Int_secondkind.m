function G = OpMat_Int_secondkind(m)
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