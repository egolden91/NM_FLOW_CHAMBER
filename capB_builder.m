%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: N - grid size
% Output: Bj matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Bout] = capB_builder(N)
B0 = zeros(N-1,N-2);
B0 = B0+[zeros(1,N-2);eye(N-2)*2];
for j=1:N-3
%     B{j} = [zeros(N-j-2);[1,zeros(1,N-j-3)]];
    B{j} = [zeros(1,N-j-2);eye(N-j-2)];
end

Bout = [B0,B];
end   