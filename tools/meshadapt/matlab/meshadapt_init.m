function meshadapt_init(ndim,smesh)
%MESHADAPT_INIT Initialisation of mesh adaptation
%
%   MESHADAPT_INIT(NDIM,SMESH)
%
%   initialises the Featflow2 mesh adaptation back-end in NDIM dimensions
%   based on the initial (coarse) mesh given by the mesh file SMESH. The
%   external data structures are accessible via the libpointer LP_MESHADAPT
%   which is stored in global memory.
%
% Input arguments:
%
%    NDIM     spatial dimension
%    SMESH    file name (including directory) of the initial mesh
%
% Author: M. Moller, TU Delft, 2014.

global lp_meshadapt

% Initialise the adaptation structure
lp_meshadapt = libpointer('t_meshAdapt');
calllib('meshadapt', 'madapt_alloc', lp_meshadapt);
calllib('meshadapt', 'madapt_init', lp_meshadapt, ...
    ndim, libpointer('stringPtr',smesh), length(smesh));
end