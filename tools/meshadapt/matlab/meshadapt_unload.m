function meshadapt_unload
%MESHADAPT_UNLOAD Unload the mesh adaptation library
%
%   MESHADAPT_UNLOAD
%
% Author: M. Moller, TU Delft, 2014.

% Unload shared library
% WARNING: Newer Matlab versions allow unloading of a shared library
%          only once. That is, loading the library, unloading it and
%          loading it again will cause a segmentation fault.
unloadlibrary meshadapt
