% Run script for solving Darcy's system (fine scale) with mixed 
% finite element method (RT0 + P0)

% clear all;

%% Set up the parameter

nx = 30;
ny = nx;

hx = 1/(nx-1);
hy = 1/(ny-1);

% load coef.mat

p = @(x)cos(pi*x(:,1)).*cos(pi*x(:,2));
v1 = @(x)-pi.*sin(pi*x(:,1)).*cos(pi*x(:,2));
v2 = @(x)-pi.*cos(pi*x(:,1)).*sin(pi*x(:,2));
% f = @(x)-2*pi^2*cos(pi*x(:,1)).*cos(pi*x(:,2));
% f = @(x,varargin)ones(size(x,1),1).*(x(:,2)>(1-2*hy)).*(x(:,1)<2*hx)-ones(size(x,1),1).*(x(:,2)<2*hy).*(x(:,1)>(1-2*hx));

%% Mesh and topology

Mesh = TProd_Mesh(0:hx:1,0:hy:1);
nElements = size(Mesh.Elements,1);
nCoordinates = size(Mesh.Coordinates,1);


% Numbering of local vertices and edges
%    4 ___ 3          ___
%     |   |          | 4 |
%     |___|        1 |___| 2
%    1     2           3

el_node = sparse([1:nElements 1:nElements 1:nElements 1:nElements],Mesh.Elements(:),[ones(1,nElements) 2*ones(1,nElements) 3*ones(1,nElements) 4*ones(1,nElements)]);
node_node = el_node'*el_node;
node_node = node_node-diag(diag(node_node));
node_node(node_node==3)=0;
node_node(node_node==8)=0;
node_node = tril(node_node);
[In,Jn,~]=find(node_node);
nEdges = numel(In);
Mesh.Vert2Edge = sparse(In,Jn,1:numel(In),nCoordinates,nCoordinates);
Mesh.Vert2Edge = Mesh.Vert2Edge + Mesh.Vert2Edge'; 
EdgeLoc = [Mesh.Vert2Edge(Mesh.Elements(:,4)+(Mesh.Elements(:,1)-1)*nCoordinates) Mesh.Vert2Edge(Mesh.Elements(:,2)+(Mesh.Elements(:,3)-1)*nCoordinates) ...
    Mesh.Vert2Edge(Mesh.Elements(:,1)+(Mesh.Elements(:,2)-1)*nCoordinates) Mesh.Vert2Edge(Mesh.Elements(:,4)+(Mesh.Elements(:,3)-1)*nCoordinates)]';

clear In Jn node_node

% plot_Mesh(Mesh,'tas');

%% Assemble system matrices of the fine saddle point problem

QuadRule_1D = gauleg(0,1,2);
QuadRule = TProd(QuadRule_1D);

% Assemble the discrete divergence operator 
% The Raviart Thomas basis function is defined such that v \cdot n = 1/|e| 
% on one edge and zero on other edges

tmp = [-ones(1,nElements);ones(1,nElements);-ones(1,nElements);ones(1,nElements)];
Ab = tmp(:); 
tmp = [1:nElements; 1:nElements; 1:nElements; 1:nElements]+nEdges;
Ib = tmp(:);
Jb = EdgeLoc(:);
EdgeLoc = EdgeLoc';

% Assemble the mass matrix (weighted by the inverse of permeability field)


low_coef = 1e-5;
coef = ones(nElements,1)*1;

margin = 0;

% horizontal channels
coef(ny/2+1+margin*(ny-1):ny-1:nElements-margin*(ny-1)) = low_coef;
coef(ny/2+margin*(ny-1):ny-1:nElements-margin*(ny-1)) = low_coef;
coef(ny/2-1+margin*(ny-1):ny-1:nElements-margin*(ny-1)) = low_coef;

coef(ny/2+4+margin*(ny-1):ny-1:nElements-margin*(ny-1)) = low_coef;
coef(ny/2+5+margin*(ny-1):ny-1:nElements-margin*(ny-1)) = low_coef;
coef(ny/2+6+margin*(ny-1):ny-1:nElements-margin*(ny-1)) = low_coef;

coef(ny/2-4+margin*(ny-1):ny-1:nElements-margin*(ny-1)) = low_coef;
coef(ny/2-5+margin*(ny-1):ny-1:nElements-margin*(ny-1)) = low_coef;
coef(ny/2-6+margin*(ny-1):ny-1:nElements-margin*(ny-1)) = low_coef;

% vertical channels
% coef((nx/2-2)*(ny-1)+(4:ny-4)) = low_coef;
% coef((nx/2-1)*(ny-1)+(4:ny-4)) = low_coef;
% coef((nx/2)*(ny-1)+(4:ny-4)) = low_coef;
% 
% coef((nx/2-5)*(ny-1)+(4:ny-4)) = low_coef;
% coef((nx/2-6)*(ny-1)+(4:ny-4)) = low_coef;
% coef((nx/2-7)*(ny-1)+(4:ny-4)) = low_coef;
% 
% coef((nx/2+3)*(ny-1)+(4:ny-4)) = low_coef;
% coef((nx/2+4)*(ny-1)+(4:ny-4)) = low_coef;
% coef((nx/2+5)*(ny-1)+(4:ny-4)) = low_coef;


Am = [(hx/hy*0.5)./coef; (hx/hy*0.5)./coef; (hy/hx*0.5)./coef; (hy/hx*0.5)./coef];

Im = [EdgeLoc(:,1); EdgeLoc(:,2); EdgeLoc(:,3); EdgeLoc(:,4)];
Jm = [EdgeLoc(:,1); EdgeLoc(:,2); EdgeLoc(:,3); EdgeLoc(:,4)];

nPorforations = 19;
A_size = nEdges+nElements+nPorforations+1;
A = sparse([Im; Ib; Jb],[Jm; Jb; Ib], [Am; Ab; Ab], A_size, A_size);

% well equation
I_v_w = nEdges+nElements+(1:nPorforations);                % porforation
I_p_w = ones(1,nPorforations)*(nEdges+nElements+nPorforations+1);  % well
I_p_wc = nEdges+(ny-1)*(nx/2-1)+ny/2+ ...
         (-(nPorforations-1)/2:(nPorforations-1)/2);  % well cell
WI = 2 * pi * coef(I_p_wc-nEdges)' * hx /log(0.28/sqrt(2)/0.05);    % well index (coef)

% 1./WI

Abw = ones(1,nPorforations);
Aw = sparse([I_v_w I_p_wc I_p_w I_v_w I_v_w], ...
            [I_v_w I_v_w I_v_w I_p_wc I_p_w], ...
            [1./WI Abw -1*Abw Abw -1*Abw], A_size, A_size);

A = A + Aw;

% A(1:nEdges, 1:nEdges)

B = A(nEdges+(1:nElements), 1:nEdges);
% Schur = B*(A(1:nEdges, 1:nEdges)\B')
% Schur(abs(Schur) < 1e-8) = 0


% 1./WI

clear Im Ib Jm Jb Am Ab;

%% Assemble right hand side

% F = zeros(nElements,1);
% nPt = numel(QuadRule.w);
% ScaledQuadRule = QuadRule.x*[hx 0;0 hy];
% for i = 1:nPt
%     F = F + QuadRule.w(i)*f([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]);
% end
% F = [zeros(nEdges,1);F*hx*hy];

F = ones(nElements,1)/nElements*(-1);
F = [zeros(nEdges,1);F;zeros(nPorforations,1);1];


%% Solve the linear system

VelocityBndDof = [EdgeLoc(1:ny-1,1); EdgeLoc(end-ny+2:end,2); EdgeLoc(1:ny-1:(ny-1)*(nx-1),3); EdgeLoc(ny-1:ny-1:(ny-1)*(nx-1),4)];
ActiveDof = setdiff([1:nEdges+nElements-1, nEdges+nElements+(1:nPorforations+1)],VelocityBndDof);

Uh=F*0;
Uh(ActiveDof) = A(ActiveDof,ActiveDof)\F(ActiveDof);
p_indices = [nEdges+1:nEdges+nElements nEdges+nElements+nPorforations+1];
Uh(p_indices) = Uh(p_indices)-sum(Uh(p_indices))*ones(nElements+1,1)/(nElements+1);

max(-Uh(p_indices))
min(-Uh(p_indices))
mean(-Uh(p_indices))

disp(['Well pressure: ', num2str(-Uh(nEdges+nElements+nPorforations+1))])
% Uh(nEdges+nElements+(1:nPorforations))
Uh = Uh(1:nEdges+nElements);
% Uh(I_p_wc)

%% Compute error
% 
% errL2v = 0;
% errL2p = 0;
% for i = 1:nPt
%     errL2v = errL2v + QuadRule.w(i)*sum((v1([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]) ...
%         -Uh(EdgeLoc(:,1))*(1/hy-QuadRule.x(i,1)/hy)-Uh(EdgeLoc(:,2))*QuadRule.x(i,1)/hy).^2);
%     errL2v = errL2v + QuadRule.w(i)*sum((v2([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]) ...
%         -Uh(EdgeLoc(:,3))*(1/hx-QuadRule.x(i,2)/hx)-Uh(EdgeLoc(:,4))*QuadRule.x(i,2)/hx).^2);
%     errL2p = errL2p + QuadRule.w(i)*sum((p([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]) ...
%         -Uh(nEdges+1:end)).^2);
% end
% errL2v = sqrt(errL2v*hx*hy);
% errL2p = sqrt(errL2p*hx*hy);
% fprintf('L2 error of the velocity is %f\n',errL2v);
% fprintf('L2 error of the pressure is %f\n',errL2p);

% Plot each component of the velocity and pressure

% v1h = [Uh(EdgeLoc(:,1));Uh(EdgeLoc(:,2));Uh(EdgeLoc(:,2));Uh(EdgeLoc(:,1))]/hy;
% v2h = [Uh(EdgeLoc(:,3));Uh(EdgeLoc(:,3));Uh(EdgeLoc(:,4));Uh(EdgeLoc(:,4))]/hx;
 
% plot_DGBFE(v1h,Mesh);
% plot_DGBFE(v2h,Mesh);
% plot_BFE(Uh(nEdges+1:end),Mesh);

% plot_Mesh(Mesh, 't')

% Numbering of local vertices and edges
%    4 ___ 3          ___
%     |   |          | 4 |
%     |___|        1 |___| 2
%    1     2           3

% plot vector field
edge_mid_pts = [
    (Mesh.Coordinates(Mesh.Elements(:,1),:)+Mesh.Coordinates(Mesh.Elements(:,4),:))/2; ...
    (Mesh.Coordinates(Mesh.Elements(:,2),:)+Mesh.Coordinates(Mesh.Elements(:,3),:))/2; ...
    (Mesh.Coordinates(Mesh.Elements(:,1),:)+Mesh.Coordinates(Mesh.Elements(:,2),:))/2; ...
    (Mesh.Coordinates(Mesh.Elements(:,4),:)+Mesh.Coordinates(Mesh.Elements(:,3),:))/2];

% Uh(:) = 1.0;
edge_v_val = [
    Uh(EdgeLoc(:,1))/hy Uh(EdgeLoc(:,3))*0; ...
    Uh(EdgeLoc(:,2))/hy Uh(EdgeLoc(:,3))*0; ...
    Uh(EdgeLoc(:,1))*0 Uh(EdgeLoc(:,3))/hx; ...
    Uh(EdgeLoc(:,1))*0 Uh(EdgeLoc(:,4))/hx];


% fig_all = figure;

% fig1 = plot_BFE(log(coef),Mesh);
% title('Log(permeability)')
% ax1 = gca;
% colorbar;
% ax1_copy = copyobj(ax1,fig_all);
% subplot(2,2,1,ax1_copy)

fig2 = plot_BFE(Uh(nEdges+1:end)*(-1),Mesh);
title('Pressure')
ax2 = gca;
colorbar;
colormap(ax2, jet);
% ax2_copy = copyobj(ax2,fig_all);
% subplot(2,2,3,ax2_copy)

fig3 = plot_Mesh(Mesh,'a');       

hold on;
quiver(edge_mid_pts(:,1), edge_mid_pts(:,2), ...
       edge_v_val(:,1), edge_v_val(:,2), 'linewidth', 1);
title('Flux')
hold off;
ax3 = gca;
% colorbar;
% ax3_copy = copyobj(ax3,fig_all);
% subplot(2,2,4,ax3_copy)

% fig4 = plot_WellCells(Mesh,I_p_wc-nEdges,'a');       
ax4 = gca;
% colorbar;
% ax4_copy = copyobj(ax4,fig_all);
% subplot(2,2,2,ax4_copy)

% 
% close(fig1)
% close(fig2)
% close(fig3)
% close(fig4)


