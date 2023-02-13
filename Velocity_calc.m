%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputs:
%       Psi - stream function values 
%       h   - step size
%   Outputs:
%       u1 - velocity in x
%       v1 - velocity in y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[u1,v1] = Velocity_calc(Psi,h)
u = [];
v = [];
% velocity calc is forward diff ((X(2)-X(1))/h)
% u calc
for qq = 1:size(Psi,2)
    u(:,qq) = diff(Psi(:,qq))/h;
end
% v calc
for qq = 1:size(Psi,1)
%     v(qq,:) = diff(Psi(qq,:))/h; %diff is the other way in y direction because psi is inverted to look convenient
    v(qq,:) = -diff(Psi(qq,:))/h; %diff is the other way in y direction because psi is inverted to look convenient
end
u1 = [u;zeros(1,length(u))];
v1 = [v,zeros(size(v,1),1)];

end