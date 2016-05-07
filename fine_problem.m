% Run script for solving Darcy's system (fine scale) with finite volume
% method in the mixed form, two point flux approximation (TPFA)

% clear all;
% tic;

%% Set up the parameter

nx = 101;
ny = nx;

hx = 1/(nx-1);
hy = 1/(ny-1);

% load coef.mat

p = @(x)cos(pi*x(:,1)).*cos(pi*x(:,2));
v1 = @(x)-pi.*sin(pi*x(:,1)).*cos(pi*x(:,2));
v2 = @(x)-pi.*cos(pi*x(:,1)).*sin(pi*x(:,2));
f = @(x)-2*pi^2*cos(pi*x(:,1)).*cos(pi*x(:,2));
% f = @(x,varargin)-ones(size(x,1),1).*(x(:,2)>(1-hy)).*(x(:,1)<hx)+ones(size(x,1),1).*(x(:,2)<hy).*(x(:,1)>(1-hx));

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

%% Assemble system matrices

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

% Assemble the mass matrix (weighted by the inverse of permeability field
coef = ones(nElements,1);
% coef = copy_spe10(Mesh);
Am = [hx/hy./coef/2; hx/hy./coef/2; hy/hx./coef/2; hy/hx./coef/2];

Im = [EdgeLoc(:,1); EdgeLoc(:,2); EdgeLoc(:,3); EdgeLoc(:,4)];
Jm = [EdgeLoc(:,1); EdgeLoc(:,2); EdgeLoc(:,3); EdgeLoc(:,4)];

% Enforce average free constraint
% Ic = (nEdges+nElements+1)*ones(nElements,1);
% Jc = (1:nElements)'+nEdges;
% Ac = ones(nElements,1);
% A = sparse([Im; Ib; Jb; Ic; Jc],[Jm; Jb; Ib; Jc; Ic], [Am; Ab; Ab; Ac; Ac]);

A = sparse([Im; Ib; Jb],[Jm; Jb; Ib], [Am; -Ab; -Ab]);
clear Im Ib Jm Jb Am Ab;

%% Assemble right hand side

F = zeros(nElements,1);
nPt = numel(QuadRule.w);
ScaledQuadRule = QuadRule.x*[hx 0;0 hy];
for i = 1:nPt
    F = F + QuadRule.w(i)*f([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]);
end
F = [zeros(nEdges,1);F*hx*hy];

%% Solve the linear system

VelocityBndDof = [EdgeLoc(1:ny-1,1); EdgeLoc(end-ny+2:end,2); EdgeLoc(1:ny-1:(ny-1)*(nx-1),3); EdgeLoc(ny-1:ny-1:(ny-1)*(nx-1),4)];
ActiveDof = setdiff(1:nEdges+nElements-1,VelocityBndDof);

Uh=F*0;
tic;
Uh(ActiveDof) = A(ActiveDof,ActiveDof)\F(ActiveDof);toc;
Uh(nEdges+1:end) = Uh(nEdges+1:end)-sum(Uh(nEdges+1:end))*ones(nElements,1)/nElements;

%% Compute error

errL2v = 0;
errL2p = 0;
for i = 1:nPt
    errL2v = errL2v + QuadRule.w(i)*sum((v1([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]) ...
        -Uh(EdgeLoc(:,1))*(1/hy-QuadRule.x(i,1)/hy)-Uh(EdgeLoc(:,2))*QuadRule.x(i,1)/hy).^2);
    errL2v = errL2v + QuadRule.w(i)*sum((v2([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]) ...
        -Uh(EdgeLoc(:,3))*(1/hx-QuadRule.x(i,2)/hx)-Uh(EdgeLoc(:,4))*QuadRule.x(i,2)/hx).^2);
    errL2p = errL2p + QuadRule.w(i)*sum((p([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]) ...
        -Uh(nEdges+1:end)).^2);
end
errL2v = sqrt(errL2v*hx*hy);
errL2p = sqrt(errL2p*hx*hy);
fprintf('L2 error of the velocity is %f\n',errL2v);
fprintf('L2 error of the pressure is %f\n',errL2p);

% toc;
% Plot each component of the velocity and pressure

% v1h = [Uh(EdgeLoc(:,1));Uh(EdgeLoc(:,2));Uh(EdgeLoc(:,2));Uh(EdgeLoc(:,1))]/hy;
% v2h = [Uh(EdgeLoc(:,3));Uh(EdgeLoc(:,3));Uh(EdgeLoc(:,4));Uh(EdgeLoc(:,4))]/hx;
% 
% plot_DGBFE(v1h,Mesh);
% plot_DGBFE(v2h,Mesh);
plot_BFE(Uh(nEdges+1:end),Mesh);
% contourf(linspace(hx/2,1,nx-1),linspace(hy/2,1,ny-1),reshape(Uh(nEdges+1:end),ny-1,nx-1),21);
% axis equal

return
%% Streamline plot

Coordinates = (Mesh.Coordinates(Mesh.Elements(:,1),:)+Mesh.Coordinates(Mesh.Elements(:,4),:)+ ...
    Mesh.Coordinates(Mesh.Elements(:,2),:)+Mesh.Coordinates(Mesh.Elements(:,3),:))/4;
U1 = (Uh(EdgeLoc(:,1))+Uh(EdgeLoc(:,2)))/2/hx;
U2 = (Uh(EdgeLoc(:,3))+Uh(EdgeLoc(:,4)))/2/hy;

space = 1;
tmp = 1:space:ny-1;
tmp = tmp(ones((nx-1)/space,1),:)';
tmp2 = (0:space:nx-2)*(ny-1);
tmp2 = tmp2(ones((ny-1)/space,1),:);
el_idx = tmp+tmp2;

Coordinates1 = el_idx*0;
Coordinates2 = el_idx*0;
U12 = el_idx*0;
U22 = el_idx*0;

for i = 1:size(Coordinates2,1)
    for j = 1:size(Coordinates2,2)
        Coordinates1(i,j) = Coordinates(el_idx(i,j),1);
        Coordinates2(i,j) = Coordinates(el_idx(i,j),2);
        U12(i,j) = U1(el_idx(i,j));
        U22(i,j) = U2(el_idx(i,j));
%         if sqrt((U12(i,j)^2+U22(i,j)^2))<.1
%             U12(i,j) = 0;
%             U22(i,j) = 0;
%         end
    end
end
% fig = quiver(Coordinates(el_idx,1),Coordinates(el_idx,2),U1(el_idx(:)),U2(el_idx(:)),'b-');
streamslice(Coordinates1,Coordinates2,U12,U22);
XMin = min(Mesh.Coordinates(:,1));
XMax = max(Mesh.Coordinates(:,1));
YMin = min(Mesh.Coordinates(:,2));
YMax = max(Mesh.Coordinates(:,2));
XLim = [XMin XMax] + .05*(XMax-XMin)*[-1 1];
YLim = [YMin YMax] + .05*(YMax-YMin)*[-1 1];
set(gca,'XLim',XLim,'YLim',YLim,'DataAspectRatio',[1 1 1]);
% plot_BFE(coef,Mesh);
plot_BFE( sqrt((U1.^2+U2.^2)),Mesh);

