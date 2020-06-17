function G = OpMat_Int_firstkind(m)
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