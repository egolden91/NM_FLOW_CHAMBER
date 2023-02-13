%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: N - grid size
%        h - step size
%        Omega - Vortex magnitute 
% Output: b vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b] = b_builder(Amat,N,h,Omega)
b = zeros(length(Amat),1);
b(ceil(N/2)) = Omega*h^2;
end