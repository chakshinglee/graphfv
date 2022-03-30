function varargout = plot_DGBFE(U,Mesh)
% PLOT_DGLFE Plot finite element solution.
%
%   PLOT_DGLFE(U,MESH) generates a plot of the finite element solution U on
%   the mesh MESH.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh.
%
%   H = PLOT_DGLFE(U,MESH) also returns the handle to the figure.
%
%   Example:
%
%   plot_DGLFE(U,MESH);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants

OFFSET = 0.05;
nElements = size(Mesh.Elements,1);

% Compute axes limits

XMin = min(Mesh.Coordinates(:,1));
XMax = max(Mesh.Coordinates(:,1));
YMin = min(Mesh.Coordinates(:,2));
YMax = max(Mesh.Coordinates(:,2));
XLim = [XMin XMax] + OFFSET*(XMax-XMin)*[-1 1];
YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];

% Generate auxiliary mesh

Coordinates = [Mesh.Coordinates(Mesh.Elements(:,1),:); Mesh.Coordinates(Mesh.Elements(:,2),:); ...
    Mesh.Coordinates(Mesh.Elements(:,3),:); Mesh.Coordinates(Mesh.Elements(:,4),:)];
Elements = [1:nElements; nElements+(1:nElements); 2*nElements+(1:nElements); 3*nElements+(1:nElements)]';

% Generate figure

if(isreal(U))
    
    % Compute color axes limits
    
    CMin = min(U);
    CMax = max(U);
    if(CMin < CMax)
        CLim = [CMin CMax] + OFFSET*(CMax-CMin)*[-1 1];
    else
        CLim = [1-OFFSET 1+OFFSET]*CMin;
    end
    
    % Plot real finite element solution

    fig = figure('Name','Discontinuous bilinear finite elements');
    patch('faces', Elements, ...
        'vertices', [Coordinates(:,1) Coordinates(:,2) U], ...
        'CData', U, ...
        'facecolor', 'interp', ...
        'edgecolor', 'none');
    set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
    set(gcf, 'Renderer','zbuffer')
    
    if(nargout > 0)
        varargout{1} = fig;
    end
    
else
    
    % Compute color axes limits
    
    CMin = min([real(U); imag(U)]);
    CMax = max([real(U); imag(U)]);
    CLim = [CMin CMax] + OFFSET*(CMax-CMin)*[-1 1];
    
    % Plot imaginary finite element solution
    
    fig_1 = figure('Name','Crouzeix-Raviart finite elements');
    patch('faces', Elements, ...
        'vertices', [Coordinates(:,1) Coordinates(:,2) real(U)], ...
        'CData', real(U), ...
        'facecolor', 'interp', ...
        'edgecolor', 'none');
    set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
    fig_2 = figure('Name','Linear finite elements');
    patch('faces', Elements, ...
        'vertices', [Coordinates(:,1) Coordinates(:,2) imag(U)], ...
        'CData', imag(U), ...
        'facecolor', 'interp', ...
        'edgecolor', 'none');
    set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
    
    if(nargout > 0)
        varargout{1} = fig_1;
        varargout{2} = fig_2;
    end
    
end

return