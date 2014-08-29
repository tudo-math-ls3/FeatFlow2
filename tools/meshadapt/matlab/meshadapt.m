function meshadapt(ndim,smesh,nref,refmax,reftol,crstol)
%MESHADAPT mesh adaptation wrapper
%
%   MESHADAPT prototypical implementation of a wrapper function which
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
    
    % Initialisation
    meshadapt_init(ndim,smesh);
        
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
    [nel,nvt] = meshadapt_data()
    
    % Get mesh from adaptation structure
    [coords,vertices,neighbours] = meshadapt_mesh()
            
    % Finalisation
    meshadapt_done();
end