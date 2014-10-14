function [nel,nvt,ndim,nnve] = meshadapt_data()
%MESHADAPT_DATA Mesh data
%
%   [NEL,NVT,NDIM] = MESHADAPT_DATA()
%
%   Returns mesh data:
%   NEL  number of elements
%   NVT  number of vertices
%   NDIM number of spatial dimensions
%
% Author: M. Moller, TU Delft, 2014.

global lp_meshadapt

% Get data from adaptation structure
nel = calllib('meshadapt', 'madapt_getnel', lp_meshadapt);
nvt = calllib('meshadapt', 'madapt_getnvt', lp_meshadapt);
ndim = calllib('meshadapt', 'madapt_getndim', lp_meshadapt);
nnve = calllib('meshadapt', 'madapt_getnnve', lp_meshadapt);