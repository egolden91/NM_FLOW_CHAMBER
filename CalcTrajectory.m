%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Arguments:
%       A: nxn matrix.
%       b: n dimensional vector.
%       omega: relaxation factor.
%       initial_guess: An initial solution guess for the solver to start with.
%       convergence_criteria: The maximum discrepancy acceptable to regard the current solution as fitting.
%   Returns:
%       Traj: solution vector of dimension n.

function [Traj,Upart,Vpart] = CalcTrajectory(PartVel,uflow,vflow,gridSize,InitLoc,dt,itrMax,plotFlag)
%% Constants 
U0    = 5e-2;      %[m/sec];
rho_p = 703;       %[kg/m^3];
D     = 100e-6;    %[m];
rho_g = 1;         %[kg/m^3];
mu    = 288.4e-7;  %[Nsec/m^2]
C1 = ((-3/4)*mu*24)/(rho_p*D^2);
C2 = D*rho_g/mu;
%%
h = gridSize;
u = uflow*U0;
v = vflow*U0;
%% Parameters 
N = 1/h;
% dt = 1/120;
% itrMax = 100000;
%%
idxrow0 = InitLoc(2);  % Y location
idxcol0 = InitLoc(1);  % X location

% idxrow0 = N/2 + 1;
% idxcol0 = N/2 + N/2 - 3;

Upart0 = PartVel(1)*max(max(u,[],'all'),max(v,[],'all'));
Vpart0 = PartVel(2)*max(max(u,[],'all'),max(v,[],'all'));
Uflow0 = u(idxrow0,idxcol0);
Vflow0 = v(idxrow0,idxcol0);

itr = 1;
while(itr<=itrMax)
    Re = C2 * sqrt((Upart0-Uflow0)^2+(Vpart0-Vflow0)^2);
    Upart(itr) = Upart0+C1*dt*(Re^0.354) * (Upart0 - Uflow0);
    Vpart(itr) = Vpart0+C1*dt*(Re^0.354) * (Vpart0 - Vflow0);

    idxrow(itr) = idxrow0+Vpart(itr) * dt/h;
    idxcol(itr) = idxcol0+Upart(itr) * dt * 2/h;

    if idxrow(itr) <= 0 || idxrow(itr) >N+1|| idxcol(itr) <= 0 || idxcol(itr) >N*2+1
        disp(['Out of grid on itr: ',num2str(itr)])
        break
    end

    idxrow0 = idxrow(itr);
    idxcol0 = idxcol(itr);

    Upart0 = Upart(itr);
    Vpart0 = Vpart(itr);
    Uflow0 = u(round(idxrow0),round(idxcol0));
    Vflow0 = v(round(idxrow0),round(idxcol0));
    itr = itr + 1;
    disp(itr)
end
xVec = -1/sqrt(2):h/sqrt(2):1/sqrt(2);
% yVec = 1/sqrt(2):-h/sqrt(2):0;
yVec = 0:h/sqrt(2):1/sqrt(2);
Traj = [xVec(round(idxcol(1:end-1)))',yVec(round(idxrow(1:end-1)))'];

if plotFlag
    figure(2)
    patch([-1/sqrt(2),1/sqrt(2),0],[0,0,1/sqrt(2)],'black','FaceAlpha',0)
    hold on
    % patch([-1/sqrt(2),1/sqrt(2),0],[0,0,-1/sqrt(2)],'black','FaceAlpha',0)
    plot(xVec(round(idxcol(1:end-1))),yVec(round(idxrow(1:end-1))),'--o','Color',[0,0.5,0],'MarkerSize',2,'MarkerFaceColor',[0,0.5,0])
    text(xVec(round(idxcol(1))),yVec(round(idxrow(1))),'S')
    text(xVec(round(idxcol(end-1))),yVec(round(idxrow(end-1))),'E')
    plot(xVec(round(idxcol(1))),yVec(round(idxrow(1))),'o','Color',[0.5,0,0],'MarkerSize',9)
    plot(xVec(round(idxcol(end-1))),yVec(round(idxrow(end-1))),'o','Color',[0.5,0,0],'MarkerSize',9)
    grid minor 
    axis equal
end
end






