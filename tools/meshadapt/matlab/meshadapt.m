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

% Perform mesh adaptation
for iref=1:nref
    
    % Get data from adaptation structure
    [nel,nvt] = meshadapt_data();
    
    % Define indicator array
    ind = zeros(nel,1); ind(1) = 2;
    
    % Perform one step of mesh adaptation
    meshadapt_step(ind,refmax,reftol,crstol);
end

% Get data from adaptation structure
[nel,nvt] = meshadapt_data();

% Get mesh from adaptation structure
[coords,vertices,neighbours] = meshadapt_mesh();

% Finalisation
meshadapt_done();
end