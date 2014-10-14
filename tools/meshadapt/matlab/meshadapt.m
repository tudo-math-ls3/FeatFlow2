function meshadapt(ndim,smesh,nref,refmax,reftol,crstol,varargin)
%MESHADAPT demo mesh adaptation wrapper
%
%   MESHADAPT(NDIM,SMESH,NREF,REFMAX,REFTOL,CRSTOL)
%   is a prototypical implementation of a wrapper function which
%   calls the mesh adaptation procedure implemented in Featflow2
%
% Input arguments:
%
%    NDIM   spatial dimension
%    SMESH  file name (including directory) of the initial mesh
%    NREF   number of refinement steps
%    REFMAX maximum number of refinement steps
%    REFTOL tolerance for refinement
%    CRSTOL tolerance for coarsening
%
% Author: M. Moller, TU Delft, 2014.

global lp_meshadapt

% Handle variable argument lists
vargin = varargin;
nargin = length(vargin);
names  = vargin(1:2:nargin);
values = vargin(2:2:nargin);

% Set default parameters
compiler = 'gcc';

% Overwrite default parameters by user-defined data
if nargin > 0
    for iargin=1:2:nargin
        switch names{iargin}
            case {'compiler'}
                compiler=values{iargin};
            otherwise
                error(['Unsupported parameter: ' names{iargin}])
        end
    end
end

% Initialisation
meshadapt_init(ndim,smesh,compiler);

% Get data from adaptation structure
[nel,nvt] = meshadapt_data();

% Get mesh from adaptation structure
[coords,vertices] = meshadapt_mesh();

% Visualise the original mesh
X = coords(1,:)';
Y = coords(2,:)';
edges = [vertices(1:2,:) vertices(2:3,:) vertices([3,1],:)];
L = sqrt(diff(X(edges),1,1).^2+diff(Y(edges),1,1).^2);
W = mean(L);

% Plot the original mesh
figure(1); clf
subplot(1,3,1)
triplot(vertices',coords(1,:),coords(2,:),'k');
axis equal tight
title('original mesh')

% Perform mesh refinement
for iref=1:nref
    
    % Create circular indicator function
    Xc = mean(X(vertices),1);
    Yc = mean(Y(vertices),1);
    phi = 0.2-sqrt((Xc-0.5).^2+(Yc-0.5).^2); % circle
    
    % Define indicator array
    ind = ones(nel,1);
    ind(abs(phi)<=W) = 2;
        
    % Perform one step of mesh adaptation
    meshadapt_step(ind,refmax,reftol,crstol);
    
    % Get data from adaptation structure
    [nel,nvt] = meshadapt_data();
    
    % Get mesh from adaptation structure
    [coords,vertices] = meshadapt_mesh();
    
    X = coords(1,:)';
    Y = coords(2,:)';
    
    W = W/2;
end

% Get data from adaptation structure
[nel,nvt] = meshadapt_data();

% Get mesh from adaptation structure
[coords,vertices] = meshadapt_mesh();

% Plot refined mesh
figure(1)
subplot(1,3,2)
triplot(vertices',coords(1,:),coords(2,:),'k');axis equal tight
title('refined mesh')

% Perform mesh re-coarsening
for iref=1:nref
    % Define indicator array
    ind = 0.1*ones(nel,1); % fully coarsen every triangle
    
    % Perform one step of mesh adaptation
    meshadapt_step(ind,refmax,reftol,crstol);
    
    % Get data from adaptation structure
    [nel,nvt] = meshadapt_data();
    
    % Get mesh from adaptation structure
    [coords,vertices] = meshadapt_mesh();
end

% plot coarsened mesh
figure(1)
subplot(1,3,3)
triplot(vertices',coords(1,:),coords(2,:),'k');axis equal tight
title('coarsened mesh')

% Finalisation
meshadapt_done();
end