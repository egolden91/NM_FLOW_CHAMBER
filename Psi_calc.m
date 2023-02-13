%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputs:
%       Amat - A matrix in Au=b
%       b    - b vector in Au=b
%       N - grid size
%       guesssor - init guess for sor solver
%       omegasor - relexation param for sor solver
%       convergence_criteria - conv creteria for sor solver
%   Outputs:
%       Psi - Stream function values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Psi] = Psi_calc(Amat,b,N,guesssor,omegasor,convergence_criteria)
%sor solver
[u2] =  sor_solver(Amat, b, omegasor, guesssor, convergence_criteria);
% solution decomposition
Z_right = zeros(N+1,N+1);
k=1;
for zz=1:N
    psi = u2(k:k+N-zz-1);
    Psi{zz} = [0;psi;0];
    k = k+N-zz;
    Z_right(:,zz) = [zeros(zz-1,1);Psi{zz}];
end
Z_left = fliplr(Z_right);
Z_left(:,end) = [];
Z = [Z_left,Z_right];
Psi = flipud(Z);

end