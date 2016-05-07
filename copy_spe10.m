function kappa = copy_spe10(Mesh)

load Udata.mat
K0=1e15*KU(:,:,:,1);%#ok

% Initialize constants

nElements = size(Mesh.Elements,1);     % Number of elements

% Generate the value of permeability at each vertex nu

kappa = zeros(nElements,1);
for i =1:nElements
    vertices = sum(Mesh.Coordinates(Mesh.Elements(i,:),:))/4;
    x = 60 + 1 - ceil(vertices(1)/20);
%     x = ceil(vertices(1)*60);
%     x = x + (x==0);
    y = ceil(vertices(2)/10);
%     y = ceil(vertices(2)*220);
%     y = y + (y==0);
    kappa(i) = K0(1,x,y);
end