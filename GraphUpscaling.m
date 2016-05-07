% Run script for solving Darcy's system (multiscale) with mixed 
% finite volume method (TPFA) for the fine scale approximation
% Modified version of Panayot's method, more favourable for parallelism 
% Use SVD to screen out the linear independency in the trace
% Extend the trace so that the basis has piecewise constant divegence

% clear all;
tic;

%% Set up the parameter

BasisPerAE = 3;
BasisPerEdge = 1*BasisPerAE;

% Load a graph from data
cd forPanayot_opte_data/
load_opte;
cd ../

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
loc_dof_size = nxAE*nyAE;
size_para = loc_dof_size*BasisPerAE;

IPp = zeros(nElements*BasisPerAE,1);   
JPp = IPp;
MPp = IPp;

size_para_b = ((nxAE-1)*nyAE+(nyAE-1)*nxAE)*(BasisPerAE-1);

IPb = zeros((Nx-1)*(Ny-1)*size_para_b,1);   
JPb = IPb;
MPb = IPb;
check = 0;
for i = 1:(Nx-1)*(Ny-1)
    el_loc = find(AE_el(i,:));
    velocity_int_dof_loc = find(AE_int_edge(i,:));
    Mloc = A(velocity_int_dof_loc,velocity_int_dof_loc);
    Bloc = A(el_loc+nEdges,velocity_int_dof_loc);
    MinvBloc = Mloc\Bloc';
    Aloc = Bloc*MinvBloc;
    Wloc = eye(numel(el_loc))*hx*hy;
    [V,lam_i]=eig(full(Aloc),full(Wloc));
    [lam_i,I]=sort(diag(lam_i),'ascend');
    V = V(:,I(1:BasisPerAE));
%     [V,~]=eigs(Aloc+Wloc,Wloc,BasisPerAE,'SM');
% check = (check+1/lam_i(BasisPerAE+1));
    MPp((i-1)*size_para+1:i*size_para) = V(:);
    tmp = el_loc(ones(1,BasisPerAE),:)';
    IPp((i-1)*size_para+1:i*size_para) = tmp(:);
    tmp = ((i-1)*BasisPerAE+1:i*BasisPerAE);
    tmp = tmp(ones(1,loc_dof_size),:);
    JPp((i-1)*size_para+1:i*size_para) = tmp(:);
    
    Bubble = MinvBloc*Wloc*V(:,2:end); 
    MPb((i-1)*size_para_b+1:i*size_para_b) = Bubble(:);
    tmp = velocity_int_dof_loc(ones(1,BasisPerAE-1),:)';
    IPb((i-1)*size_para_b+1:i*size_para_b) = tmp(:);
    tmp = ((i-1)*(BasisPerAE-1)+1:i*(BasisPerAE-1));
    tmp = tmp(ones(1,numel(velocity_int_dof_loc)),:);
    JPb((i-1)*size_para_b+1:i*size_para_b) = tmp(:);
        
end
toc;
%% Constructing spectral basis for the velocity 

tic;
loc_dof_size = ((nxAE-1)*nyAE+(nyAE-1)*nxAE)*2+nyAE;
half_loc_dof_size = (nxAE-1)*nyAE+(nyAE-1)*nxAE+nyAE;
Ioff = zeros((Nx*(Ny-1)+Ny*(Nx-1))*loc_dof_size*BasisPerEdge,1);   
Joff = Ioff;
Moff = Ioff;
size_para_v = loc_dof_size*(2*BasisPerAE-1);

cnt = 0;
for i = 1:Nx*(Ny-1)+(Nx-1)*Ny
    [AE_no,~,~] = find(CEdgeLoc==i);
    if numel(AE_no) == 2
        % Need to separate the two coarse blocks so as to impose the average conditions and rhs easier 
        snap_dof_loc = intersect(find(AE_bnd_edge(AE_no(1),:)),find(AE_bnd_edge(AE_no(2),:)));
        trace = zeros(numel(snap_dof_loc),2*(BasisPerAE-1));
        for k = 1:2
            el_ext_loc = find(AE_el_ext(AE_no(k),:));
            velocity_ext_int_dof_loc = find(AE_int_edge_ext(AE_no(k),:));
            Mloc = A(velocity_ext_int_dof_loc,velocity_ext_int_dof_loc);
            Bloc = A(el_ext_loc+nEdges,velocity_ext_int_dof_loc);
%             if (AE_no(k)==25) 
%                 hi=100;
%             end
            MinvBloc = Mloc\Bloc';
            Aloc = Bloc*MinvBloc;
            Wloc = eye(numel(el_ext_loc))*hx*hy;
            [V,lam_i]=eig(full(Aloc),full(Wloc));
            [lam_i,I]=sort(diag(lam_i),'ascend');% assert(abs(lam_i(1))<1e-10,'first eigenvalue is not zero');
            V = V(:,I(1:BasisPerAE));
            check = max(check,1/lam_i(BasisPerAE+1));
            %     [V,~]=eigs(Aloc+Wloc,Wloc,BasisPerAE,'SM');
            snap_marker = ismember(velocity_ext_int_dof_loc,snap_dof_loc);
            assert(norm(velocity_ext_int_dof_loc(snap_marker)-snap_dof_loc)==0,'snap_marker is not correct');
            trace(:,(k-1)*(BasisPerAE-1)+(1:BasisPerAE-1)) = MinvBloc(snap_marker,:)*Wloc*V(:,2:end);
        end
        el_loc_1 = find(AE_el(AE_no(1),:));
        velocity_int_dof_loc_1 = find(AE_int_edge(AE_no(1),:));
        
        el_loc_2 = find(AE_el(AE_no(2),:));
        velocity_int_dof_loc_2 = find(AE_int_edge(AE_no(2),:));
        
        velocity_int_dof_loc = [velocity_int_dof_loc_1 velocity_int_dof_loc_2];
        velocity_dof_loc = [snap_dof_loc velocity_int_dof_loc]; 
        pressure_dof_loc = nEdges+[el_loc_1 el_loc_2];
        int_dof_loc = [velocity_dof_loc pressure_dof_loc];
        
        size_v = numel(velocity_dof_loc);
        Aloc = A(int_dof_loc,int_dof_loc);
        Floc = zeros(numel(int_dof_loc),1);
        Floc(size_v+1:end,1) = [ones(nxAE*nyAE,1);-ones(nxAE*nyAE,1)]*hx*hy;
        Floc = Floc*hx*hy;
        
        ActiveDofLoc = 1:size_v+2*nxAE*nyAE-1;
        Uloc_int = Aloc(ActiveDofLoc,ActiveDofLoc)\Floc(ActiveDofLoc);
%         
%         const_one = ones(numel(snap_dof_loc),1);
        const_one = Uloc_int(1:numel(snap_dof_loc));
        for ti = 1:size(trace,2)
            trace_i = trace(:,ti);
            trace(:,ti) = trace_i-(trace_i'*const_one)/(const_one'*const_one)*const_one;
            if norm(trace(:,ti))<1e-10
                trace(:,ti) = 0;
            end
        end
%         trace2 = [ones(numel(snap_dof_loc),1) trace]; 
%         trace2 = [Uloc_int(1:numel(snap_dof_loc)) trace]; 
        trace2 = trace;
        [U,S,V] = svd(trace2);
        S = diag(S);
        S = S/max(S);
%         trace = trace2*V(:,S>1e-8);
        trace = U(:,S>1e-8);
%         trace = [ones(numel(snap_dof_loc),1) trace]; %#ok
        trace = [Uloc_int(1:numel(snap_dof_loc)) trace]; %#ok
        if sum(S>1e-8)>(BasisPerEdge-1)
            trace = trace(:,1:BasisPerEdge);
        end
%         trace = ones(numel(snap_dof_loc),1);
%         trace = [trace1/2 trace2];
        
        velocity_int_dof_loc = [find(AE_int_edge(AE_no(1),:)) find(AE_int_edge(AE_no(2),:))];
        velocity_dof_loc = [velocity_int_dof_loc snap_dof_loc]; 
        pressure_dof_loc = nEdges+[find(AE_el(AE_no(1),:)) find(AE_el(AE_no(2),:))];
        size_v = numel(velocity_int_dof_loc);
        int_dof_loc = [velocity_int_dof_loc pressure_dof_loc];
        Aloc = A(int_dof_loc,int_dof_loc);
        Floc = -A(int_dof_loc,snap_dof_loc)*trace;
        Floc2 = [zeros(size_v,size(Floc,2)); ones(nxAE*nyAE,1)*sum(Floc(size_v+1:size_v+nxAE*nyAE,:))/(nxAE*nyAE); ...
            ones(nxAE*nyAE,1)*sum(Floc(size_v+nxAE*nyAE+1:size_v+2*nxAE*nyAE,:))/(nxAE*nyAE)];
        Floc3 = Floc - Floc2;
        
        ActiveDofLoc = [1:size_v size_v+2:size_v+2*nxAE*nyAE-1];
        Uloc_int = Aloc(ActiveDofLoc,ActiveDofLoc)\Floc3(ActiveDofLoc,:);
        Uloc = [Uloc_int(1:size_v,:);trace];
        
        size_para_loc = numel(Uloc);
        Moff((i-1)*size_para_v+(1:size_para_loc)) = Uloc(:);
        tmp = velocity_dof_loc(ones(1,size(Uloc,2)),:)';
        Ioff((i-1)*size_para_v+(1:size_para_loc)) = tmp(:);
        tmp = cnt+(1:size(Uloc,2));
        tmp = tmp(ones(1,loc_dof_size),:);
        Joff((i-1)*size_para_v+(1:size_para_loc)) = tmp(:);   
        cnt = cnt+size(Uloc,2);
    end
end
Joff = Joff(Ioff>0);
Moff = Moff(Ioff>0);
Ioff = Ioff(Ioff>0);
P_off = sparse([Ioff; IPb; IPp+nEdges], [Joff+(Nx-1)*(Ny-1)*(2*BasisPerAE-1); JPb+(Nx-1)*(Ny-1)*BasisPerAE; ... 
    JPp], [Moff; MPb; MPp]);

toc;

%% Solve the coarse scale linear system

M = A(1:nEdges,1:nEdges);
B = A(nEdges+1:end,1:nEdges);

Pp = P_off(nEdges+1:end,1:(Nx-1)*(Ny-1)*BasisPerAE);
Pu = P_off(1:nEdges,(Nx-1)*(Ny-1)*BasisPerAE+1:end);

% A = [speye(nEdges) M\B'; B zeros(nElements)];
% A = [M B'; B zeros(nElements)];
% A_off2 = Pp'*(B*Pu*((Pu'*M*Pu)\Pu'*B'))*Pp;
A_off2 = Pp'*B*(M\B')*Pp;
F_off = P_off'*F;
Uh_off2 = -Pp(:,2:end)*(A_off2(2:end,2:end)\F_off(2:(Nx-1)*(Ny-1)*BasisPerAE));
Uh_off2 = Uh_off2-sum(Uh_off2)*ones(nElements,1)/nElements;
Uh_off3 = [-M\(B'*Uh_off2);Uh_off2];

% tic;
A_off = P_off'*A*P_off;
F_off = P_off'*F;

% activeDof = [];
% for i = 1:BasisPerEdge
%     activeDof = union(activeDof,[(CEdgeLoc(1:Ny-1,1)-1)*BasisPerEdge+i (CEdgeLoc(end-(Ny-1)+1:end,2)-1)*BasisPerEdge+i ...
%         (CEdgeLoc(1:Ny-1:nCElements,3)-1)*BasisPerEdge+i (CEdgeLoc(Ny-1:Ny-1:nCElements,4)-1)*BasisPerEdge+i]);
% end
% activeDof = setdiff(1:size(A_off,2)-1,activeDof);
% activeDof = setdiff(1:size(A_off,2),[activeDof; size(A_off,2)-BasisPerAE+1]);
activeDof = 2:size(A_off,2);

% A_off(abs(A_off)<1e-14)=0;
Uh_off = A_off(activeDof,activeDof)\F_off(activeDof); 
% Uh_off = minres(A_off(activeDof,activeDof),F_off(activeDof),1e-9,100000);
Uh_off = P_off(:,activeDof)*Uh_off;
Uh_off(nEdges+1:end) = Uh_off(nEdges+1:end)-sum(Uh_off(nEdges+1:end))*ones(nElements,1)/nElements;
Uh_off = Uh_off3;
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

plot_DGBFE(v1h,Mesh);
plot_DGBFE(v2h,Mesh);
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