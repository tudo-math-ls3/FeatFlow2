function meshadapt_done
%MESHADAPT_DONE Finalisation of mesh adaptation
%
%   MESHADAPT_DONE
%
%   deallocates all externally allocated memory and stored as libpointer
%   LP_MESHADAPT in global memory and finalises the Featflow2 back-end.
%
% Author: M. Moller, TU Delft, 2014.

global lp_meshadapt

% Finalisation
calllib('meshadapt', 'madapt_done', lp_meshadapt);
calllib('meshadapt', 'madapt_dealloc', lp_meshadapt);
clear global lp_meshadapt

% Clean up the storage management, finish
calllib('meshadapt', 'storage_done', [])

% Unload shared library
% WARNING: Newer Matlab versions allow unloading of a shared library
%          only once. That is, loading the library, unloading it and
%          loading it again will cause a segmentation fault.
unloadlibrary meshadapt
