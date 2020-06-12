function Obj = PathObj_Fourier_Flat(Z,m,j,Jval)
    
    % Unpack Fourier Coefficients
    alp1 = Z(1:m,:);
    bet1 = Z(m+1:2*m,:);
    alp2 = Z(2*m+1:3*m,:);
    bet2 = Z(3*m+1:4*m,:);
    alp3 = Z(4*m+1:5*m,:);
    bet3 = Z(5*m+1:6*m,:);
    
    
    % Flat Frame Objective
    if Jval == 1 
        Obj = sum(j'.*(abs(alp1)+abs(bet1)))^2 + sum(j'.*(abs(alp2)+abs(bet2)))^2 + sum(j'.*(abs(alp3)+abs(bet3)))^2;
    elseif Jval == 2
        Obj = sum(j'.^2.*(abs(alp1)+abs(bet1)))^2 + sum(j'.^2.*(abs(alp2)+abs(bet2)))^2 + sum(j'.^2.*(abs(alp3)+abs(bet3)))^2;
    else
        error("Unknown Objective Function Selected, Continuing with Minimum Velocity");
        Obj = sum(j'.*(abs(alp1)+abs(bet1)))^2 + sum(j'.*(abs(alp2)+abs(bet2)))^2 + sum(j'.*(abs(alp3)+abs(bet3)))^2;
    end
end 