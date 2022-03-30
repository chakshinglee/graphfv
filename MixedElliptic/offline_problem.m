% Run script for solving Darcy's system with mixed
% generalized multiscale finite element method

clear all;
tic;

%% Set up the parameter

Nx = 11;
Ny = Nx;

nx = 161;
ny = nx;

BasisPerEdge = 1;

Hx = 1/(Nx-1);
Hy = 1/(Ny-1);

hx = 1/(nx-1);
hy = 1/(ny-1);

nxAE = (nx-1)/(Nx-1);
nyAE = (ny-1)/(Ny-1);

load coef.mat
coef = ones((nx-1)*(ny-1),1);
p = @(x)cos(pi*x(:,1)).*cos(pi*x(:,2));
v1 = @(x)-pi.*sin(pi*x(:,1)).*cos(pi*x(:,2));
v2 = @(x)-pi.*cos(pi*x(:,1)).*sin(pi*x(:,2));
% f = @(x)-2*pi^2*cos(pi*x(:,1)).*cos(pi*x(:,2));
f = @(x,varargin)ones(size(x,1),1).*(x(:,2)>(1-Hy)).*(x(:,1)<Hx)-ones(size(x,1),1).*(x(:,2)<Hy).*(x(:,1)>(1-Hx));

%% Fine mesh and topology

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
    Mesh.Vert2Edge(Mesh.Elements(:,1)+(Mesh.Elements(:,2)-1)*nCoordinates) Mesh.Vert2Edge(Mesh.Elements(:,4)+(Mesh.Elements(:,3)-1)*nCoordinates)];

% plot_Mesh(Mesh,'tas');

CMesh = TProd_Mesh(0:Hx:1,0:Hy:1);
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

%% Agglomerated Elements (coarse mesh) and topology

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

% coef = ones(nElements,1);
Am = [sum((1/hy-QuadRule.x(:,1)/hy).*(1/hy-QuadRule.x(:,1)/hy).*QuadRule.w)./coef; ...
    sum((QuadRule.x(:,1)/hy).*(QuadRule.x(:,1)/hy).*QuadRule.w)./coef; ...
    sum((1/hx-QuadRule.x(:,2)/hx).*(1/hx-QuadRule.x(:,2)/hx).*QuadRule.w)./coef; ...
    sum((QuadRule.x(:,2)/hx).*(QuadRule.x(:,2)/hx).*QuadRule.w)./coef; ...
    sum((1/hy-QuadRule.x(:,1)/hy).*(QuadRule.x(:,1)/hy).*QuadRule.w)./coef; ...
    sum((1/hy-QuadRule.x(:,1)/hy).*(QuadRule.x(:,1)/hy).*QuadRule.w)./coef; ...
    sum((1/hx-QuadRule.x(:,2)/hx).*(QuadRule.x(:,2)/hx).*QuadRule.w)./coef; ...
    sum((1/hx-QuadRule.x(:,2)/hx).*(QuadRule.x(:,2)/hx).*QuadRule.w)./coef]*hx*hy;

Im = [EdgeLoc(:,1); EdgeLoc(:,2); EdgeLoc(:,3); EdgeLoc(:,4); EdgeLoc(:,1); EdgeLoc(:,2); EdgeLoc(:,3); EdgeLoc(:,4)];
Jm = [EdgeLoc(:,1); EdgeLoc(:,2); EdgeLoc(:,3); EdgeLoc(:,4); EdgeLoc(:,2); EdgeLoc(:,1); EdgeLoc(:,4); EdgeLoc(:,3)];

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

VelocityBndDof = [EdgeLoc(1:ny-1,1) EdgeLoc(end-ny+2:end,2) EdgeLoc(1:ny-1:(ny-1)*(nx-1),3) EdgeLoc(ny-1:ny-1:(ny-1)*(nx-1),4)];
ActiveDof = setdiff(1:nEdges,VelocityBndDof(:));

B = A(nEdges+1:end,ActiveDof);
M = A(ActiveDof,ActiveDof);
Lapl = B*(M\B');
Lapl = (Lapl+Lapl')/2;
eigs(Lapl,2,'SM')
return;


Uh=F*0;
Uh(ActiveDof) = A(ActiveDof,ActiveDof)\F(ActiveDof);
Uh(nEdges+1:end) = Uh(nEdges+1:end)-sum(Uh(nEdges+1:end))*ones(nElements,1)/nElements;

%% Compute the norms of the fine solution

normL2v = sqrt(Uh(1:nEdges)'*A(1:nEdges,1:nEdges)*Uh(1:nEdges));
normL2p = norm(Uh(nEdges+1:end))*sqrt(hx*hy);


%% Constructing offline space basis for the veloccity space

loc_dof_size = ((nxAE-1)*nyAE+(nyAE-1)*nxAE)*2+nyAE;
half_loc_dof_size = (nxAE-1)*nyAE+(nyAE-1)*nxAE+nyAE;
Ioff = zeros((Nx*(Ny-1)+Ny*(Nx-1))*loc_dof_size*BasisPerEdge,1);
Joff = Ioff;
Moff = Ioff;

coefInverseSum = zeros(nEdges,1);
coefInverseSum(EdgeLoc(:,1)) = coefInverseSum(EdgeLoc(:,1)) + 1./coef;
coefInverseSum(EdgeLoc(:,2)) = coefInverseSum(EdgeLoc(:,2)) + 1./coef;
coefInverseSum(EdgeLoc(:,3)) = coefInverseSum(EdgeLoc(:,3)) + 1./coef;
coefInverseSum(EdgeLoc(:,4)) = coefInverseSum(EdgeLoc(:,4)) + 1./coef;

% loop over each coarse edge

for i = 1:Nx*(Ny-1)+(Nx-1)*Ny
    [AE_no,~,~] = find(CEdgeLoc==i);
    if numel(AE_no) == 2
        % Extract the local dofs in each coarse neighborhood (union of two coarse elements) of an interior coarse edge
        % Extract separately the two coarse elements so as to impose the average conditions and rhs easier
        
        el_loc_1 = find(AE_el(AE_no(1),:));
        velocity_int_dof_loc_1 = find(AE_int_edge(AE_no(1),:));
        
        el_loc_2 = find(AE_el(AE_no(2),:));
        velocity_int_dof_loc_2 = find(AE_int_edge(AE_no(2),:));
        
        velocity_int_dof_loc = [velocity_int_dof_loc_1 velocity_int_dof_loc_2];
        pressure_dof_loc = nEdges+[el_loc_1 el_loc_2];
        int_dof_loc = [velocity_int_dof_loc pressure_dof_loc];
        snap_dof_loc = intersect(find(AE_bnd_edge(AE_no(1),:)),find(AE_bnd_edge(AE_no(2),:)));
        
        % Extract local (saddle point) matrix and assemble local rhs (with compatibility condition satisfied)
        
        size_para = numel(velocity_int_dof_loc);
        Aloc = A(int_dof_loc,int_dof_loc);
        Floc = -A(int_dof_loc,snap_dof_loc);
        Floc2 = [zeros(size_para,size(Floc,2)); ones(nxAE*nyAE,1)*sum(Floc(size_para+1:size_para+nxAE*nyAE,:))/(nxAE*nyAE); ...
            ones(nxAE*nyAE,1)*sum(Floc(size_para+nxAE*nyAE+1:size_para+2*nxAE*nyAE,:))/(nxAE*nyAE)];
        Floc3 = Floc - Floc2;
        
        % Solve for the snapshot basis functions
        
        ActiveDofLoc = [1:size_para size_para+2:size_para+2*nxAE*nyAE-1];
        Uloc_int = Aloc(ActiveDofLoc,ActiveDofLoc)\Floc3(ActiveDofLoc,:);
        Uloc = [Uloc_int(1:size_para,:);speye(numel(snap_dof_loc))];
        
        % Solve generalized eigenvalue problem on the snapshot space to produce offline basis
        
        Aloc_snap = Uloc'*A([velocity_int_dof_loc snap_dof_loc],[velocity_int_dof_loc snap_dof_loc])*Uloc ...
            +2*ones(numel(snap_dof_loc))/(nxAE*nyAE);
        Aloc_snap = (Aloc_snap+Aloc_snap')/2;
        Xloc_snap = diag(coefInverseSum(snap_dof_loc))/nxAE;
        [V,lam_i]=eig(full(Aloc_snap),full(Xloc_snap));
        [lam_i,I]=sort(diag(lam_i),'descend');
        V = V(:,I(1:BasisPerEdge));
        
        Uloc = Uloc*V;
        
        % Store the indices and entries of the offline basis
        
        size_para = loc_dof_size*BasisPerEdge;
        Moff((i-1)*size_para+1:i*size_para) = Uloc(:);
        tmp = [velocity_int_dof_loc(ones(1,BasisPerEdge),:)'; snap_dof_loc(ones(1,BasisPerEdge),:)'];
        Ioff((i-1)*size_para+1:i*size_para) = tmp(:);
        tmp = ((i-1)*BasisPerEdge+1:i*BasisPerEdge);
        tmp = tmp(ones(1,loc_dof_size),:);
        Joff((i-1)*size_para+1:i*size_para) = tmp(:);
    else
        % Extract the local dofs in the coarse element attached to a coarse edge on the domain boudary
        
        el_loc = find(AE_el(AE_no,:));
        velocity_int_dof_loc = find(AE_int_edge(AE_no,:));
        
        pressure_dof_loc = nEdges+el_loc;
        int_dof_loc = [velocity_int_dof_loc pressure_dof_loc];
        snap_dof_loc = intersect(find(AE_bnd_edge(AE_no,:)),VelocityBndDof(:,CEdgeLoc(AE_no,:)==i))';
        
        % Extract local (saddle point) matrix and assemble local rhs (with compatibility condition satisfied)
        
        size_para = numel(velocity_int_dof_loc);
        Aloc = A(int_dof_loc,int_dof_loc);
        Floc = -A(int_dof_loc,snap_dof_loc);
        Floc2 = [zeros(size_para,size(Floc,2)); ones(nxAE*nyAE,1)*sum(Floc(size_para+1:end,:))/(nxAE*nyAE)];
        Floc3 = Floc - Floc2;
        
        % Solve for the snapshot basis functions
        
        Uloc_int = Aloc(1:end-1,1:end-1)\Floc3(1:end-1,:);
        Uloc = [Uloc_int(1:size_para,:);speye(numel(snap_dof_loc))];
        
        % Solve generalized eigenvalue problem on the snapshot space to produce offline basis
        
        Aloc_snap = Uloc'*A([velocity_int_dof_loc snap_dof_loc],[velocity_int_dof_loc snap_dof_loc])*Uloc ...
            +ones(numel(snap_dof_loc))/(nxAE*nyAE);
        Aloc_snap = (Aloc_snap+Aloc_snap')/2;
        Xloc_snap = diag(coefInverseSum(snap_dof_loc))/nxAE;
        [V,lam_i]=eig(full(Aloc_snap),full(Xloc_snap));
        [lam_i,I]=sort(diag(lam_i),'descend');
        V = V(:,I(1:BasisPerEdge));
        
        Uloc = Uloc*V;
        
        % Store the indices and entries of the offline basis
        
        size_para = loc_dof_size*BasisPerEdge;
        Moff((i-1)*size_para+(1:half_loc_dof_size*BasisPerEdge)) = Uloc(:);
        tmp = [velocity_int_dof_loc(ones(1,BasisPerEdge),:)'; snap_dof_loc(ones(1,BasisPerEdge),:)'];
        Ioff((i-1)*size_para+(1:half_loc_dof_size*BasisPerEdge)) = tmp(:);
        tmp = ((i-1)*BasisPerEdge+1:i*BasisPerEdge);
        tmp = tmp(ones(1,half_loc_dof_size),:);
        Joff((i-1)*size_para+(1:half_loc_dof_size*BasisPerEdge)) = tmp(:);
    end
end
Joff = Joff(Ioff>0);
Moff = Moff(Ioff>0);
Ioff = Ioff(Ioff>0);

% Since coarse pressure basis is piecewise constant on the coarse grid, we can simply use
% the agglomerated element to fine element map to define the coarse pressure basis

[Ioff2, Joff2, Moff2] = find(AE_el);

% Construct the coarse-to-fine mapping for both the velocity and pressure spaces

P_off = sparse([Ioff; Joff2+nEdges], [Joff; Ioff2+(Nx*(Ny-1)+(Nx-1)*Ny)*BasisPerEdge], [Moff; Moff2]);

%% Assemble and solve the coarse scale linear system

A_off = P_off'*A*P_off;
F_off = P_off'*F;

activeDof = [];
for i = 1:BasisPerEdge
    activeDof = union(activeDof,[(CEdgeLoc(1:Ny-1,1)-1)*BasisPerEdge+i (CEdgeLoc(end-(Ny-1)+1:end,2)-1)*BasisPerEdge+i ...
        (CEdgeLoc(1:Ny-1:nCElements,3)-1)*BasisPerEdge+i (CEdgeLoc(Ny-1:Ny-1:nCElements,4)-1)*BasisPerEdge+i]);
end
activeDof = setdiff(1:size(A_off,2)-1,activeDof);

tic;
Uh_off = A_off(activeDof,activeDof)\F_off(activeDof);
Uh_off = P_off(:,activeDof)*Uh_off;
Uh_off(nEdges+1:end) = Uh_off(nEdges+1:end)-sum(Uh_off(nEdges+1:end))*ones(nElements,1)/nElements;
toc;
Uh = Uh-Uh_off;

%% Compute error

errL2v = sqrt(Uh(1:nEdges)'*A(1:nEdges,1:nEdges)*Uh(1:nEdges))/normL2v;
errL2p = norm(Uh(nEdges+1:end))*sqrt(hx*hy)/normL2p;

fprintf('Relative L2 error of the velocity is %f\n',errL2v);
fprintf('Relative L2 error of the pressure is %f\n',errL2p);

% Plot each component of the velocity and pressure

Uh = Uh_off;
v1h = [Uh(EdgeLoc(:,1));Uh(EdgeLoc(:,2));Uh(EdgeLoc(:,2));Uh(EdgeLoc(:,1))]/hy;
v2h = [Uh(EdgeLoc(:,3));Uh(EdgeLoc(:,3));Uh(EdgeLoc(:,4));Uh(EdgeLoc(:,4))]/hx;

plot_DGBFE(v1h,Mesh);
plot_DGBFE(v2h,Mesh);
plot_BFE(Uh(nEdges+1:end),Mesh);


