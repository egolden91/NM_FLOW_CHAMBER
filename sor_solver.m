%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputs:
%       A - nxn matrix.
%       b - n dimensional vector.
%       omega -  relaxation factor.
%       initial_guess - An initial solution guess for the solver to start with.
%       convergence_criteria -  The maximum discrepancy acceptable to regard the current solution as fitting.
%   Outputs:
%       x- solution vector of dimension n.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[x] =  sor_solver(A, b, omega, initial_guess, convergence_criteria)
b = b(:);
step = 0;
x = initial_guess(:);
residual = norm((A*x) - b);  % Initial residual
while residual > convergence_criteria
    for i = 1:size(A,1)
        sigma = 0;
        for j = 1: size(A,2)
            if j ~= i
                sigma=sigma+A(i,j)*x(j);
            end
        end
        x(i) = (1 - omega) * x(i) + (omega/A(i, i))*(b(i) - sigma);
    end
    residual = norm(A*x - b);
    step = step + 1;
end
end