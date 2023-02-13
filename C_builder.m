%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: N - grid size
% Output: Cj matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cout] = C_builder(N)
for j=1:N-3
%     C{j} = [[zeros(1,N-j-1),1];zeros(N-j-2,N-j)];
    C{j} = [zeros(N-j-1,1),eye(N-j-1)];
end
C{N-2} = [0,1];
Cout = C;
end