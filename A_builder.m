%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: N - grid size
% Output: Aj matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Aout] = A_builder(N)
A0 = zeros(N-1,N-1);
A0 = A0+eye(N-1)*(-4);
for ii= 1:N-2
    A0(ii,ii+1) = 1;
    A0(ii+1,ii) = 1;
end

for j = 1:N-3
    A{j} = zeros(N-j-1)+eye(N-j-1)*(-4);
    for ii= 1:N-j-2
        A{j}(ii,ii+1) = 1;
        A{j}(ii+1,ii) = 1;
    end
end
A{N-2} = -4;
Aout = [A0,A];
end