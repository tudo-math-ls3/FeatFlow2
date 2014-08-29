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

    [nel,nvt,ndim] = meshadapt_data()
    
    % Get coordinates
    p_coords = libpointer('doublePtr', zeros(ndim,nvt));
    calllib('meshadapt', 'meshadaptbase_MOD_madapt_getvertexcoords', ...
            lp_meshadapt, p_coords);
    coords = p_coords.Value;

    % Get vertices
    p_vertices = libpointer('int32Ptr', zeros(4,nel));
    calllib('meshadapt', 'meshadaptbase_MOD_madapt_getverticesatelement', ...
            lp_meshadapt, p_vertices);
    vertices = p_vertices.Value;
    
    % Get adjacencies
    p_neighbours = libpointer('int32Ptr', zeros(4,nel));
    calllib('meshadapt', 'meshadaptbase_MOD_madapt_getneighboursatelement', ...
            lp_meshadapt, p_neighbours);
    neighbours = p_neighbours.Value;