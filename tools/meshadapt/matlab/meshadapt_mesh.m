function [coords,vertices,neighbours] = meshadapt_mesh()
%MESHADAPT_MESH Mesh
%
%   [COORDS,VERTICES,NEIGHBOURS] = MESHADAPT_DATA()
%
%   Returns mesh data:
%   COORDS     vertex coordinates
%   VERTICES   vertices-at-element
%   NEIGHBOURS elements-at-element
%
% Author: M. Moller, TU Delft, 2014.

global lp_meshadapt

[nel,nvt,ndim,nnve] = meshadapt_data();

if nargout==0, return, end

% Get coordinates
p_coords = libpointer('doublePtr', zeros(ndim,nvt));
calllib('meshadapt', 'madapt_getvertexcoords', ...
    lp_meshadapt, p_coords);
coords = p_coords.Value;

if nargout==1, return, end

% Get vertices
p_vertices = libpointer('int32Ptr', zeros(nnve,nel));
calllib('meshadapt', 'madapt_getverticesatelement', ...
    lp_meshadapt, p_vertices);
vertices = p_vertices.Value;

if nargout==2, return, end

% Get adjacencies
p_neighbours = libpointer('int32Ptr', zeros(nnve,nel));
calllib('meshadapt', 'madapt_getneighboursatelement', ...
    lp_meshadapt, p_neighbours);
neighbours = p_neighbours.Value;