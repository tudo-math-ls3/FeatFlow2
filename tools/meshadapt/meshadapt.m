% Load shared library (if required)
if ~libisloaded('meshadapt'), 
    [notfound,warnings] = ...
        loadlibrary('libmeshadapt.so','libmeshadapt.h', ...
        'addheader', 'types', ...
        'alias','meshadapt');
end

% DEBUG: show library functions and structures
libfunctions meshadapt

% Initialise system-wide settings
calllib('meshadapt', 'fsystem_MOD_sys_init_simple')

% Initialise the output system
calllib('meshadapt', 'genoutput_MOD_output_init_simple')

% Initialise the FEAT 2.0 storage management
calllib('meshadapt', 'storage_MOD_storage_init', 100, 100, [])

% Initialisation
calllib('meshadapt','meshadaptbase_MOD_madapt_init', 2, 'mesh/bench1', []);

% Finalisation
calllib('meshadapt','meshadaptbase_MOD_madapt_done', []);

% Clean up the storage management, finish
calllib('meshadapt', 'storage_MOD_storage_done', [])

% Unload shared library
unloadlibrary meshadapt