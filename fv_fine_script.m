Grid.Nx=10; Grid.hx=1/Grid.Nx;
Grid.Ny=10; Grid.hy=1/Grid.Ny;
Grid.Nz=1; Grid.hz=1/Grid.Nz;
% Grid.K=ones(3,Grid.Nx,Grid.Ny);
Grid.K=exp(5*smooth3(smooth3(randn(3,Grid.Nx,Grid.Ny))));
N=Grid.Nx*Grid.Ny*Grid.Nz; q=zeros(N,1); q([1 N])=[1 -1];
tic;
P=TPFA_mixed(Grid,Grid.K,q);
toc;
contourf(linspace(Grid.hx/2,1,Grid.Nx),linspace(Grid.hy/2,1,Grid.Ny),P,11);
axis square