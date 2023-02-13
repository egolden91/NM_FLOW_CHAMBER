%% Hyper Parameters
clear all
% i runs from 0 to N-i
% j runs from 0 to N-j
tic
h = 1/64;
N = 1/h;
Omega = 800;
%%Stream Function and Velocity field Calculation
[Aj] = A_builder(N);
% B matrices
[Bj] = capB_builder(N);
% C matrices
[Cj] = C_builder(N);
% Main matrix
[Amat] = MainMat_builder(Aj,Bj,Cj,N);
% b vector
[b] = b_builder(Amat,N,h,Omega);
% Calculation of the stream function
% sor method params calc
Cjac = eye(length(Amat)) - inv(diag(diag(Amat)))*Amat;
mu = max(abs(eig(Cjac)));
OMEGA = 1+(mu/(1+sqrt(1-mu^2))^2);
guessSOR = zeros(1,length(b));
omegaSOR = OMEGA;
convSOR = 1e-8;
Psi = Psi_calc(Amat,b,N,guessSOR,omegaSOR,convSOR);
% Calculation of the velocity field
[u,v] = Velocity_calc(Psi,h);
toc
%% Trajectory calculation of particals 
PartVel = [0.2,0;4,0;10,0;40,0;50,0;100,0]*-1;
PartVel = [50,0;100,0]*-1;
% InitLoc = [2*N - N/2 ,N/2];
InitLoc = [2*N - N/2 - 1 ,N/2];
dt = 1/1000;
itrMax = 10000000; 
for zz = 1:size(PartVel,1)
    [Traj{zz},Upart{zz},Vpart{zz}] = CalcTrajectory(PartVel(zz,:),u,v,h,InitLoc,dt,itrMax,0);
end
%% PLOTS
% Plot velocity field 
figure(1)
xVec = -1/sqrt(2):h/sqrt(2):1/sqrt(2);
% yVec = 1/sqrt(2):-h/sqrt(2):0;
yVec = 0:h/sqrt(2):1/sqrt(2);
[X,Y] = meshgrid(xVec,yVec);
figure(1);
contourf(X,Y,Psi)
colorbar
hold on
% quiver(X,Y,flipud(u),flipud(v))
quiver(X,Y,u,v)
patch([-1/sqrt(2),1/sqrt(2),0],[0,0,1/sqrt(2)],'black','FaceAlpha',0)
title('Stream contour and velocity field \Omega_0 = 800')
xlabel('X')
ylabel('Y')
% patch([-1/sqrt(2),1/sqrt(2),0],[0,0,-1/sqrt(2)],'black','FaceAlpha',0)
axis equal
grid minor
%%
%Plot trajectory
for zz = 1:size(PartVel,1)
    figure(zz)
    xVec = -1/sqrt(2):h/sqrt(2):1/sqrt(2);
    % yVec = 1/sqrt(2):-h/sqrt(2):0;
    yVec = 0:h/sqrt(2):1/sqrt(2);
    [X,Y] = meshgrid(xVec,yVec);
    contourf(X,Y,Psi)
    colorbar
    hold on
    quiver(X,Y,u,v)
    patch([-1/sqrt(2),1/sqrt(2),0],[0,0,1/sqrt(2)],'black','FaceAlpha',0)
    plot(Traj{zz}(:,1),Traj{zz}(:,2),'--o','Color',[0,0.5,0],'MarkerSize',3,'MarkerFaceColor',[0,0.5,0])
    text(Traj{zz}(1,1),Traj{zz}(1,2),'S')
    text(Traj{zz}(end,1),Traj{zz}(end,2),'E')
    plot(Traj{zz}(1,1),Traj{zz}(1,2),'o','Color',[0.5,0,0],'MarkerSize',9)
    plot(Traj{zz}(end,1),Traj{zz}(end,2),'o','Color',[0.5,0,0],'MarkerSize',9)
    axis equal
    grid minor
    title(['Partical trajectory ',' M(',num2str(PartVel(zz,1)),',0)'])
    xlabel('X')
    ylabel('Y')
end
%% ANALYSIS
% %% Analysis Grid Size
% PSI = {};
% U = {};
% V = {};
% h = [1/16,1/32,1/64,1/100];
% Omega = 800;
% for i = 1:length(h)
%     N = 1/h(i);
%     % Stream Function and Velocity field Calculation
%     [Aj] = A_builder(N);
%     % B matrices
%     [Bj] = capB_builder(N);
%     % C matrices
%     [Cj] = C_builder(N);
%     % Main matrix
%     [Amat] = MainMat_builder(Aj,Bj,Cj,N);
%     [b] = b_builder(Amat,N,h(i),Omega);
%     %     condNum(i) = cond(Amat);% sor method params calc
%     Cjac = eye(length(Amat)) - inv(diag(diag(Amat)))*Amat;
%     mu = max(abs(eig(Cjac)));
%     OMEGA = 1+(mu/(1+sqrt(1-mu^2))^2);
%     guessSOR = zeros(1,length(b));
%     omegaSOR = OMEGA;
%     convSOR = 1e-8;
%     [~,PSI{i}] = Psi_calc(Amat,b,N,guessSOR,omegaSOR,convSOR);
%     [U{i},V{i}] = Velocity_calc(PSI{i},h(i));
% end
% for ii = 1:length(h)
%     xVec = -1/sqrt(2):h(ii)/sqrt(2):1/sqrt(2);
%     % yVec = 1/sqrt(2):-h/sqrt(2):0;
%     yVec = 0:h(ii)/sqrt(2):1/sqrt(2);
%     [X,Y] = meshgrid(xVec,yVec);
%     figure(4)
%     contour(X,Y,PSI{ii},'Color',[rand(1,3)])
%     grid minor
%     hold on
% end
% %% dt analysis
% Traj1 = {};
% Upart1= {};
% Vpart1= {};
% PartVel = [50,0];
% itrMax = 100000;
% InitLoc = [N/2+1 ,N/2];
% dt = [1/10 1/32 1/64 1/100 1/200 1/400 1/500 1/800 1/1000 1/2000 1/10000];
% for ii = 1:length(dt)
%    [Traj1{ii},Upart1{ii},Vpart1{ii}] = CalcTrajectory(PartVel,u,v,h,InitLoc,dt(ii),itrMax,0);
% end
% 
% for jj = 5:length(Traj1)
%     figure(11)
%     plot(Traj1{jj}(:,1),Traj1{jj}(:,2),'--','LineWidth',0.5)
%     hold on
%     grid minor
% end
% patch([-1/sqrt(2),1/sqrt(2),0],[0,0,1/sqrt(2)],'black','FaceAlpha',0)