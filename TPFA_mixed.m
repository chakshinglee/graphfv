function[P,V]=TPFA_mixed(Grid,K,q)
% Compute transmissibilities by harmonic averaging.
Nx=Grid.Nx; Ny=Grid.Ny; Nz=Grid.Nz; N=Nx*Ny*Nz;
hx=Grid.hx; hy=Grid.hy; hz=Grid.hz;
L = K.^(-1);
tx = 2*hy*hz/hx; TX = zeros(Nx+1,Ny,Nz);
ty = 2*hx*hz/hy; TY = zeros(Nx,Ny+1,Nz);
tz = 2*hx*hy/hz; TZ = zeros(Nx,Ny,Nz+1);
TX(2:Nx,:,:) = tx./(L(1,1:Nx-1,:,:)+L(1,2:Nx,:,:));
TY(:,2:Ny,:) = ty./(L (2,:,1: Ny-1,:)+L(2,:,2:Ny,:));
TZ(:,:,2:Nz) = tz./(L (3,:,:,1: Nz-1)+L(3,:,:,2:Nz));

% Assemble TPFA discretization matrix.
x1 = reshape(TX(1:Nx,:,:),N,1); x2 = reshape(TX(2:Nx+1,:,:),N,1);
y1 = reshape(TY(:,1:Ny,:),N,1); y2 = reshape(TY(:,2:Ny+1,:),N,1);
z1 = reshape(TZ(:,:,1:Nz),N,1); z2 = reshape(TZ(:,:,2:Nz+1),N,1);
DiagVecs = [-z2,-y2,-x2,x1+x2+y1+y2+z1+z2,-x1,-y1,-z1];
DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
A = spdiags(DiagVecs,DiagIndx,N,N);
A(1,1) = A(1,1)+sum(Grid.K(:,1,1,1));

% Assemble v_{ij} = t_{ij}(u_i - u_j) 
helper_x = reshape(1:(Nx+1)*Ny*Nz,Nx+1,Ny,Nz);
helper_y = reshape(1:Nx*(Ny+1)*Nz,Nx,Ny+1,Nz);
helper_z = reshape(1:Nx*Ny*(Nz+1),Nx,Ny,Nz+1);
IX1 = reshape(helper_x(1:Nx,:,:),1,N); IX2 = reshape(helper_x(2:Nx+1,:,:),1,N);
IY1 = reshape(helper_y(:,1:Ny,:),1,N); IY2 = reshape(helper_y(:,2:Ny+1,:),1,N);
IZ1 = reshape(helper_z(:,:,1:Nz),1,N); IZ2 = reshape(helper_z(:,:,2:Nz+1),1,N);

os_x = (Nx+1)*Ny*Nz; % offset for x-direction
os_y = Nx*(Ny+1)*Nz; % offset for y-direction

BT = sparse([IX1 IX2 os_x+IY1 os_x+IY2 os_x+os_y+IZ1 os_x+os_y+IZ2], ...
    [1:N 1:N 1:N 1:N 1:N 1:N],[-x1; x2; -y1; y2; -z1; z2]);  
B = sparse([1:N 1:N 1:N 1:N 1:N 1:N],[IX1 os_x+IY1 os_x+os_y+IZ1 ...
    IX2 os_x+IY2 os_x+os_y+IZ2],[-ones(1,N*3) ones(1,N*3)]);  

norm(A-B*BT,'fro')

% Solve linear system and extract interface fluxes.
u = A\q;
% u = (B*BT)\q;
P = reshape(u,Nx,Ny,Nz);
V.x = zeros(Nx+1,Ny,Nz);
V.y = zeros(Nx,Ny+1,Nz);
V.z = zeros(Nx,Ny,Nz+1);
V.x(2:Nx,:,:) = (P(1:Nx-1,:,:)-P(2:Nx,:,:)).*TX(2:Nx,:,:);
V.y(:,2:Ny,:) = (P(:,1:Ny-1,:)-P(:,2:Ny,:)).*TY(:,2:Ny,:);
V.z(:,:,2:Nz) = (P(:,:,1:Nz-1)-P(:,:,2:Nz)).*TZ(:,:,2:Nz);