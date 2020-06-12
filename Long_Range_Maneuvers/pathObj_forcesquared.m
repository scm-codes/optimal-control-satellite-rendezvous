function obj = pathObj_forcesquared(U)
% obj = pathObj_forcesquared(U)
%

% Unpack U
ux = U(1,:);
uy = U(2,:);
uz = U(3,:);

% Force Squared objective
obj = ux.^2+uy.^2 + uz.^2;
end