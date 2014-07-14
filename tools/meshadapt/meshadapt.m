% Load shared library (if required)
if ~libisloaded('meshadapt'),
    
    if ispc,
        loadlibrary('libmeshadapt.dll','meshadapt.h', ...
                    'addheader', 'types', ...
                    'alias','meshadapt');
    elseif ismac,
        loadlibrary('libmeshadapt.dylib','meshadapt.h', ...
                    'addheader', 'types', ...
                    'alias','meshadapt');
    elseif ~ismac && isunix,
        loadlibrary('libmeshadapt.so','meshadapt.h', ...
                    'addheader', 'types', ...
                    'alias','meshadapt');
    else
        error('Unsupported operating system!')
    end
end

% DEBUG: show library functions and structures
libfunctions meshadapt -full

% Initialise system-wide settings
calllib('meshadapt', 'fsystem_MOD_sys_init_simple')

% Initialise the output system
calllib('meshadapt', 'genoutput_MOD_output_init_simple')

% Initialise the FEAT 2.0 storage management
calllib('meshadapt', 'storage_MOD_storage_init', 100, 100, [])

% Initialisation
t_meshAdapt = libpointer('t_meshAdapt');
calllib('meshadapt','meshadaptbase_MOD_madapt_alloc',t_meshAdapt);
%calllib('meshadapt','meshadaptbase_MOD_madapt_init', t_meshAdapt, 2, ...
%        libpointer('stringPtr','mesh/bench1'));

% Finalisation
calllib('meshadapt','meshadaptbase_MOD_madapt_done', t_meshAdapt);
calllib('meshadapt','meshadaptbase_MOD_madapt_dealloc', t_meshAdapt);
clear t_meshAdapt

% Clean up the storage management, finish
calllib('meshadapt', 'storage_MOD_storage_done', [])

% Unload shared library
unloadlibrary meshadapt