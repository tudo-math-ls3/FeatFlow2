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
%    NDIM    spatial dimension
%    SMESH   file name (including directory) of the initial mesh
%
% Author: M. Moller, TU Delft, 2014.
    
    global lp_meshadapt
    
    % Load shared library (if required)
    if ~libisloaded('meshadapt'),
        if ispc,
            loadlibrary('libmeshadapt.dll', 'meshadapt.h', ...
                        'addheader', 'types', ...
                        'alias','meshadapt');
        elseif ismac,
            loadlibrary('libmeshadapt.dylib', 'meshadapt.h', ...
                        'addheader', 'types', ...
                        'alias','meshadapt');
        elseif ~ismac && isunix,
            loadlibrary('libmeshadapt.so', 'meshadapt.h', ...
                        'addheader', 'types', ...
                        'alias','meshadapt');
        else
            error('Unsupported operating system!')
        end
    end

    % Initialise system-wide settings
    calllib('meshadapt', 'fsystem_MOD_sys_init_simple')

    % Initialise the output system
    calllib('meshadapt', 'genoutput_MOD_output_init_simple')

    % Initialise the FEAT 2.0 storage management
    calllib('meshadapt', 'storage_MOD_storage_init', 100, 100, [])

    % Initialise the adaptation structure
    lp_meshadapt = libpointer('t_meshAdapt');
    calllib('meshadapt', 'meshadaptbase_MOD_madapt_alloc', lp_meshadapt);
    calllib('meshadapt', 'meshadaptbase_MOD_madapt_init', lp_meshadapt, ...
            ndim, libpointer('stringPtr',smesh), length(smesh));
end