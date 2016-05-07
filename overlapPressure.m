% Run script for solving Darcy's system (multiscale) with mixed 
% finite volume method (TPFA) for the fine scale approximation
% Modified version of Panayot's method, more favourable for parallelism 
% Use SVD to screen out the linear independency in the trace
% Extend the trace so that the basis has piecewise constant divegence

% clear all;
tic;

%% Set up the parameter

Nx = 7;
Ny = 23;

nx = 61;
ny = 221;

Lx = 1200;
Ly = 2200;

BasisPerAE = 80;
BasisPerEdge = 1*BasisPerAE;

Hx = Lx/(Nx-1);
Hy = Ly/(Ny-1);

hx = Lx/(nx-1);
hy = Ly/(ny-1);

nxAE = (nx-1)/(Nx-1);
nyAE = (ny-1)/(Ny-1);

% load coef.mat

p = @(x)cos(pi*x(:,1)).*cos(pi*x(:,2));
v1 = @(x)-pi.*sin(pi*x(:,1)).*cos(pi*x(:,2));
v2 = @(x)-pi.*cos(pi*x(:,1)).*sin(pi*x(:,2));
% f = @(x)-2*pi^2*cos(pi*x(:,1)).*cos(pi*x(:,2));
f = @(x,varargin)ones(size(x,1),1).*(x(:,2)>(Ly-Hy)).*(x(:,1)<Hx)-ones(size(x,1),1).*(x(:,2)<Hy).*(x(:,1)>(Lx-Hx));
% perm = @(x)100 - 99.999*(( ((x(:,1)-.75).^2+(x(:,2)-.6).^2) <.0225 ) + ( ((x(:,1)-1.25).^2+(x(:,2)-.4).^2) <.0225 ) );
% perm = @(x)100 - 99.999*(( ((x(:,1)-.95).^2+(x(:,2)-.7).^2) <.0225 ) + ( ((x(:,1)-1.05).^2+(x(:,2)-.3).^2) <.0225 ) );

%% Mesh and topology

Mesh = TProd_Mesh(0:hx:Lx,0:hy:Ly);
nElements = size(Mesh.Elements,1);
nCoordinates = size(Mesh.Coordinates,1);

% Numbering of loca vertices and edges
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
    Mesh.Vert2Edge(Mesh.Elements(:,1)+(Mesh.Elements(:,2)-1)*nCoordinates) Mesh.Vert2Edge(Mesh.Elements(:,4)+(Mesh.Elements(:,3)-1)*nCoordinates)];

% plot_Mesh(Mesh,'tas');

CMesh = TProd_Mesh(0:Hx:Lx,0:Hy:Ly);
nCElements = size(CMesh.Elements,1);
nCCoordinates = size(CMesh.Coordinates,1);

Cel_node = sparse([1:nCElements 1:nCElements 1:nCElements 1:nCElements],CMesh.Elements(:),[ones(1,nCElements) 2*ones(1,nCElements) 3*ones(1,nCElements) 4*ones(1,nCElements)]);
Cnode_node = Cel_node'*Cel_node;
Cnode_node = Cnode_node-diag(diag(Cnode_node));
Cnode_node(Cnode_node==3)=0;
Cnode_node(Cnode_node==8)=0;
Cnode_node = tril(Cnode_node);
[In,Jn,~]=find(Cnode_node);
CMesh.Vert2Edge = sparse(In,Jn,1:numel(In),nCCoordinates,nCCoordinates);
CMesh.Vert2Edge = CMesh.Vert2Edge + CMesh.Vert2Edge';

CEdgeLoc = [CMesh.Vert2Edge(CMesh.Elements(:,4)+(CMesh.Elements(:,1)-1)*nCCoordinates) CMesh.Vert2Edge(CMesh.Elements(:,2)+(CMesh.Elements(:,3)-1)*nCCoordinates) ...
    CMesh.Vert2Edge(CMesh.Elements(:,1)+(CMesh.Elements(:,2)-1)*nCCoordinates) CMesh.Vert2Edge(CMesh.Elements(:,4)+(CMesh.Elements(:,3)-1)*nCCoordinates)];
clear In Jn node_node

%% Agglomerates and topology
tic;
el_edge = sparse([1:nElements 1:nElements 1:nElements 1:nElements],EdgeLoc(:),ones(1,4*nElements));
tmp = mod((1:nElements),ny-1); tmp(tmp==0) = ny-1;
% Mesh.ElemFlag = ceil(tmp/(nyAE)) + (Ny-1)*(ceil(ceil((1:nElements)/(ny-1))/(nxAE))-1);
AE_el = sparse(1:nElements, ceil(tmp/(nyAE)) + (Ny-1)*(ceil(ceil((1:nElements)/(ny-1))/(nxAE))-1), ones(1,nElements));
AE_el = AE_el';

AE_edge = AE_el*el_edge;
AE_int_edge = AE_edge;
AE_int_edge(AE_int_edge~=2)=0;
AE_bnd_edge = AE_edge-AE_int_edge;
EdgeLoc = EdgeLoc';
clear tmp;

el_el = el_edge*el_edge';
AE_el_ext = AE_el*el_el;
AE_el_ext(AE_el_ext>0) = 1;
AE_int_edge_ext = AE_el_ext*el_edge;
AE_int_edge_ext(AE_int_edge_ext~=2)=0;

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
% coef = ones(nElements,1);
coef = copy_spe10(Mesh);
% load ../MixedElliptic/coef.mat
% coef = perm((Mesh.Coordinates(Mesh.Elements(:,1),:)+Mesh.Coordinates(Mesh.Elements(:,2),:)+...
%     Mesh.Coordinates(Mesh.Elements(:,3),:)+Mesh.Coordinates(Mesh.Elements(:,4),:))/4);

Am = [hx/hy./coef/2; hx/hy./coef/2; hy/hx./coef/2; hy/hx./coef/2];

Im = [EdgeLoc(:,1); EdgeLoc(:,2); EdgeLoc(:,3); EdgeLoc(:,4)];
Jm = [EdgeLoc(:,1); EdgeLoc(:,2); EdgeLoc(:,3); EdgeLoc(:,4)];

% Am = [sum((1/hy-QuadRule.x(:,1)/hy).*(1/hy-QuadRule.x(:,1)/hy).*QuadRule.w)./coef; ...
% sum((QuadRule.x(:,1)/hy).*(QuadRule.x(:,1)/hy).*QuadRule.w)./coef; ...
% sum((1/hx-QuadRule.x(:,2)/hx).*(1/hx-QuadRule.x(:,2)/hx).*QuadRule.w)./coef; ...
% sum((QuadRule.x(:,2)/hx).*(QuadRule.x(:,2)/hx).*QuadRule.w)./coef; ...
% sum((1/hy-QuadRule.x(:,1)/hy).*(QuadRule.x(:,1)/hy).*QuadRule.w)./coef; ...
% sum((1/hy-QuadRule.x(:,1)/hy).*(QuadRule.x(:,1)/hy).*QuadRule.w)./coef; ...
% sum((1/hx-QuadRule.x(:,2)/hx).*(QuadRule.x(:,2)/hx).*QuadRule.w)./coef; ...
% sum((1/hx-QuadRule.x(:,2)/hx).*(QuadRule.x(:,2)/hx).*QuadRule.w)./coef]*hx*hy;
% 
% Im = [EdgeLoc(:,1); EdgeLoc(:,2); EdgeLoc(:,3); EdgeLoc(:,4); EdgeLoc(:,1); EdgeLoc(:,2); EdgeLoc(:,3); EdgeLoc(:,4)];
% Jm = [EdgeLoc(:,1); EdgeLoc(:,2); EdgeLoc(:,3); EdgeLoc(:,4); EdgeLoc(:,2); EdgeLoc(:,1); EdgeLoc(:,4); EdgeLoc(:,3)];

% Enforce average free constraint
% Ic = (nEdges+nElements+1)*ones(nElements,1);
% Jc = (1:nElements)'+nEdges;
% Ac = ones(nElements,1);
% A = sparse([Im; Ib; Jb; Ic; Jc],[Jm; Jb; Ib; Jc; Ic], [Am; Ab; Ab; Ac; Ac]);

A = sparse([Im; Ib; Jb],[Jm; Jb; Ib], [Am; Ab; Ab]);
clear Im Ib Jm Jb Am Ab;

%% Assemble right hand side

F = zeros(nElements,1);
nPt = numel(QuadRule.w);
ScaledQuadRule = QuadRule.x*[hx 0;0 hy];
for i = 1:nPt
    F = F + QuadRule.w(i)*f([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]);
end
F = [zeros(nEdges,1);F*hx*hy];

%% Solve the fine scale linear system

VelocityBndDof = [EdgeLoc(1:ny-1,1); EdgeLoc(end-ny+2:end,2); EdgeLoc(1:ny-1:(ny-1)*(nx-1),3); EdgeLoc(ny-1:ny-1:(ny-1)*(nx-1),4)];
ActiveDof = setdiff(1:nEdges+nElements-1,VelocityBndDof(:));

Uh=F*0;
Uh(ActiveDof) = A(ActiveDof,ActiveDof)\F(ActiveDof);
Uh(nEdges+1:end) = Uh(nEdges+1:end)-sum(Uh(nEdges+1:end))*ones(nElements,1)/nElements;

%% Compute the norms of the fine solution

normL2v = sqrt(Uh(1:nEdges)'*A(1:nEdges,1:nEdges)*Uh(1:nEdges));
normL2p = norm(Uh(nEdges+1:end))*sqrt(hx*hy);
toc;

%% Constructing spectral basis for the pressure
tic;
loc_dof_size = (nxAE+2)*(nyAE+2);
size_para = loc_dof_size*BasisPerAE;

IPp = zeros(2*nElements*BasisPerAE,1);   
JPp = IPp;
MPp = IPp;

for i = 1:(Nx-1)*(Ny-1)
    el_loc = find(AE_el_ext(i,:));
    velocity_int_dof_loc = find(AE_int_edge_ext(i,:));
    Mloc = A(velocity_int_dof_loc,velocity_int_dof_loc);
    Bloc = A(el_loc+nEdges,velocity_int_dof_loc);
%     MinvBloc = spdiags(1./diag(Mloc),0,size(Mloc,1),size(Mloc,1))*Bloc';
    MinvBloc = Mloc\Bloc';
    Aloc = Bloc*MinvBloc;
    Wloc = eye(numel(el_loc))*hx*hy;
    [V,lam_i]=eig(full(Aloc),full(Wloc));
    [lam_i,I]=sort(diag(lam_i),'ascend');
    V = V(:,I(1:BasisPerAE));
    
    size_tmp = numel(V);
    MPp((i-1)*size_para+(1:size_tmp)) = V(:);
    tmp = el_loc(ones(1,BasisPerAE),:)';
    IPp((i-1)*size_para+(1:size_tmp)) = tmp(:);
    tmp = ((i-1)*BasisPerAE+1:i*BasisPerAE);
    tmp = tmp(ones(1,numel(el_loc)),:);
    JPp((i-1)*size_para+(1:size_tmp)) = tmp(:);
end
JPp = JPp(IPp>0);
MPp = MPp(IPp>0);
IPp = IPp(IPp>0);
Pp = sparse(IPp, JPp, MPp);

toc;

%% Solve the coarse scale linear system

M = A(1:nEdges,1:nEdges);
B = A(nEdges+1:end,1:nEdges);

% Pp = P_off(nEdges+1:end,1:(Nx-1)*(Ny-1)*BasisPerAE);
% Pu = P_off(1:nEdges,(Nx-1)*(Ny-1)*BasisPerAE+1:end);

% A = [speye(nEdges) M\B'; B zeros(nElements)];
% A = [M B'; B zeros(nElements)];
% A_off2 = Pp'*(B*Pu*((Pu'*M*Pu)\Pu'*B'))*Pp;
A_off = Pp'*B*(M\B')*Pp;
F_off = Pp'*F(nEdges+1:end);
Uh_off = -Pp(:,1:end)*(A_off(1:end,1:end)\F_off(1:end));
Uh_off = Uh_off-sum(Uh_off)*ones(nElements,1)/nElements;
Uh_off = [-M\(B'*Uh_off);Uh_off];

diff = Uh-Uh_off;

%% Compute error

errL2v = sqrt(diff(1:nEdges)'*A(1:nEdges,1:nEdges)*diff(1:nEdges));
errL2p = norm(diff(nEdges+1:end))*sqrt(hx*hy);
errL2v = errL2v/normL2v;
errL2p = errL2p/normL2p;

fprintf('Relative L2 error of the velocity is %f\n',errL2v);
fprintf('Relative L2 error of the pressure is %f\n',errL2p);

% Plot each component of the velocity and pressure

Uh = Uh_off;
% Uh = P_off(:,3);
v1h = [Uh(EdgeLoc(:,1));Uh(EdgeLoc(:,2));Uh(EdgeLoc(:,2));Uh(EdgeLoc(:,1))]/hy;
v2h = [Uh(EdgeLoc(:,3));Uh(EdgeLoc(:,3));Uh(EdgeLoc(:,4));Uh(EdgeLoc(:,4))]/hx;

% plot_DGBFE(v1h,Mesh);
% plot_DGBFE(v2h,Mesh);
plot_BFE(Uh(nEdges+1:end),Mesh);

return
%% Plot streamline

% Coordinates = [Mesh.Coordinates(Mesh.Elements(:,1),:)+Mesh.Coordinates(Mesh.Elements(:,4),:); ...
%     Mesh.Coordinates(Mesh.Elements(:,2),:)+Mesh.Coordinates(Mesh.Elements(:,3),:); ...
%     Mesh.Coordinates(Mesh.Elements(:,1),:)+Mesh.Coordinates(Mesh.Elements(:,2),:); ...
%     Mesh.Coordinates(Mesh.Elements(:,3),:)+Mesh.Coordinates(Mesh.Elements(:,4),:)]/2;
% U1 = [Uh(EdgeLoc(:,1)); Uh(EdgeLoc(:,2)); 0*Uh(EdgeLoc(:,3)); 0*Uh(EdgeLoc(:,4))];
% U2 = [0*Uh(EdgeLoc(:,1)); 0*Uh(EdgeLoc(:,2)); Uh(EdgeLoc(:,3)); Uh(EdgeLoc(:,4))];
Coordinates = (Mesh.Coordinates(Mesh.Elements(:,1),:)+Mesh.Coordinates(Mesh.Elements(:,4),:)+ ...
    Mesh.Coordinates(Mesh.Elements(:,2),:)+Mesh.Coordinates(Mesh.Elements(:,3),:))/4;
U1 = (Uh(EdgeLoc(:,1))+Uh(EdgeLoc(:,2)))/2;
U2 = (Uh(EdgeLoc(:,3))+Uh(EdgeLoc(:,4)))/2;

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
        if sqrt((U12(i,j)^2+U22(i,j)^2))<.1
            U12(i,j) = 0;
            U22(i,j) = 0;
        end
    end
end
% fig = quiver(Coordinates(el_idx,1),Coordinates(el_idx,2),U1(el_idx(:)),U2(el_idx(:)),'b-');
% streamslice(Coordinates1,Coordinates2,U12,U22);
XMin = min(Mesh.Coordinates(:,1));
XMax = max(Mesh.Coordinates(:,1));
YMin = min(Mesh.Coordinates(:,2));
YMax = max(Mesh.Coordinates(:,2));
XLim = [XMin XMax] + .05*(XMax-XMin)*[-1 1];
YLim = [YMin YMax] + .05*(YMax-YMin)*[-1 1];
% set(gca,'XLim',XLim,'YLim',YLim,'DataAspectRatio',[1 1 1]);
% plot_BFE(log10(coef),Mesh);
plot_BFE( log10(sqrt((U1.^2+U2.^2))),Mesh);