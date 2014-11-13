function meshadapt_load(compiler)
%MESHADAPT_INIT Load mesh adaptation library
%
%   MESHADAPT_LOAD(COMPILER)
%
%   load the Featflow2 mesh adaptation back-end.
%
% Input arguments:
%
%    COMPILER name of the compiler suite
%
% Author: M. Moller, TU Delft, 2014.

% Set default parameters
headerfile = ['meshadapt_' compiler];
protofile  = ['meshadapt_' compiler];
if ispc,
    libfile = 'libmeshadapt.dll';
elseif ismac,
    libfile = 'libmeshadapt.dylib';
elseif ~ismac && isunix,
    libfile = 'libmeshadapt.so';
else
    error('Unsupported operating system!');
end

% Load shared library (if required)
if ~libisloaded('meshadapt'),
    
    % Check if prototype file exists for this type of compiler
    if exist ([protofile '.m'], 'file') == 2,
        
        % Prototype file exists
        try
            % Try to load the shared library based on the prototype file
            loadlibrary(libfile, str2func(protofile),...
                'alias','meshadapt');
        catch
            % Most likely the file meshadapt_thunk_COMPUTER.LIBEXT does not
            % exist. Hence, do a dummy load of the library using the
            % original header file to create the aforementioned library.
            try
                loadlibrary(libfile, [headerfile '.h'], ...
                    'addheader', 'types', ...
                    'alias','meshadapt',...
                    'mfilename','dummy.m')
                unloadlibrary('meshadapt');
                
                % Remove the dumme prototype file
                delete('dummy.m')
                
                % Load the shared library based on the prototype file
                loadlibrary(libfile, str2func(protofile),...
                    'alias','meshadapt');
            catch
                error('Unable to load shared library')
            end
        end
        
    else
        
        % Prototype file does not exist. Load the shared library using
        % the existing header files, generate a prototypical prototype
        % file, ask the user to adjust it to the particular compiler
        % and stop continuation of the program.
        try
            loadlibrary(libfile, [headerfile '.h'], ...
                'addheader', 'types', ...
                'alias','meshadapt',...
                'mfilename',protofile);
            unloadlibrary('meshadapt');
            
            disp(['A template prototype file ' protofile ...
                ' has been created from the header file ' headerfile '.']);
            disp(['Please adjust the prototype file as follows:']);
            disp('');
            disp('1) Replace the line');
            disp('   ThunkLibName=fullfile(MfilePath,"meshadapt_thunk_XYY");');
            disp('by');
            disp('   ThunkLibName=fullfile(MfilePath,["meshadapt_thunk_" lower(computer)]);');
            disp('');
            disp('2) Add aliases for each library function.');
            
        catch
            error('An error occured while loading the shared library using the header files.')
        end
        error('Mesh adaptation cannot be performed.')
    end
end

% Initialise system-wide settings
calllib('meshadapt', 'sys_init_simple')

% Initialise the output system
calllib('meshadapt', 'output_init_simple')

% Initialise the FEAT 2.0 storage management
calllib('meshadapt', 'storage_init', 100, 100, [])