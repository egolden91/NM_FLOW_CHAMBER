%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputs:
%       Aj - Aj matrix cell
%       Bj - Bj matrix cell
%       Cj - Cj matrix cell
%       N - Grid size
%   Outputs:
%       Amat - full A matrix in Au=b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Amat] = MainMat_builder(Aj,Bj,Cj,N)
MatWidth = length(Aj)*(1+length(Aj))/2;
Amat = [];
Amat = [Aj{1},Bj{1},zeros(size(Aj{1},1),MatWidth-size(Aj{1},2)-size(Bj{1},2))];
x = 0;
for qq = 2:N-2
    Arow = zeros(size(Aj{qq},1),MatWidth);
    CABmat = horzcat(Cj{qq-1},Aj{qq},Bj{qq});
    r = size(CABmat,1);
    c = size(CABmat,2);
    if qq ~=2
        x(qq-1) = N - qq + 2;
        Arow(1:r,sum(x)+1:sum(x)+c) = CABmat;
    else
        Arow(1:r,1:c) = CABmat;
    end
    Amat = [Amat;Arow];
end
Amat = [Amat;[zeros(size(Aj{end},1),MatWidth-size(Aj{end},2)-size(Cj{end},2)),Cj{end},Aj{end}]];

end