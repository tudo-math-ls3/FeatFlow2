function meshadapt_init(ndim,smesh,compiler)
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
%    COMPILER name of the compiler suite
%
% Author: M. Moller, TU Delft, 2014.
    
    global lp_meshadapt
    
    headerfile = ['meshadapt_' compiler];
    protofile  = ['meshadapt_' compiler];
    
    % Load shared library (if required)
    if ~libisloaded('meshadapt'),
        % Check if prototype file exists
        if exist ([protofile '.m'], 'file') == 2,
            % Prototype file exists. So we can load the shared library and continue with
            % mesh adaptation
            if ispc,
                loadlibrary('libmeshadapt.dll', str2func(protofile));
            elseif ismac,
                loadlibrary('libmeshadapt.dylib', str2func(protofile));
            elseif ~ismac && isunix,
                loadlibrary('libmeshadapt.so', str2func(protofile),...
                            'alias','meshadapt');
            else
                error('Unsupported operating system!')
            end
        else
            % Prototype file does not exist. Load the shared library,
            % generate a prototypical prototype file and ask the user to
            % adjust it to the particular compiler.
            if ispc,
                loadlibrary('libmeshadapt.dll', [headerfile '.h'], ...
                            'addheader', 'types', ...
                            'alias','meshadapt',...
                            'mfilename',protofile);
            elseif ismac,
                loadlibrary('libmeshadapt.dylib', [headerfile '.h'], ...
                            'addheader', 'types', ...
                            'alias','meshadapt',...
                            'mfilename',protofile);
            elseif ~ismac && isunix,
                loadlibrary('libmeshadapt.so', [headerfile '.h'], ...
                            'addheader', 'types', ...
                            'alias','meshadapt',...
                            'mfilename',protofile);
            else
                error('Unsupported operating system!')
            end
            disp(['Created template prototype file ' protofile ...
                  ' from header file ' headerfile '.']);
            disp(['Please adjust it to your compiler configuration.']);
            unloadlibrary('meshadapt');
            return
        end
    end

    % Initialise system-wide settings
    calllib('meshadapt', 'sys_init_simple')

    % Initialise the output system
    calllib('meshadapt', 'output_init_simple')

    % Initialise the FEAT 2.0 storage management
    calllib('meshadapt', 'storage_init', 100, 100, [])

    % Initialise the adaptation structure
    lp_meshadapt = libpointer('t_meshAdapt');
    calllib('meshadapt', 'madapt_alloc', lp_meshadapt);
    calllib('meshadapt', 'madapt_init', lp_meshadapt, ...
            ndim, libpointer('stringPtr',smesh), length(smesh));
end