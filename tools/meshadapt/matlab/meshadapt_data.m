function [varargout] = meshadapt_data(useTria,varargin)
%MESHADAPT_DATA Mesh data
%
%   [VARARGOUT] = MESHADAPT_DATA(USETRIA,VARARGIN)
%
%   Input arguments:
%   USETRIA    retrieve data from triangulation structure
%   VARARGIN   list of strings defining the output
%
%   Returns mesh data:
%   NEL        number of elements
%   NVT        number of vertices
%   NDIM       number of spatial dimensions
%   NNVE       maximm number of vertices per element
%   COORDS     vertex coordinates
%   VERTICES   vertices-at-element
%   NEIGHBOURS elements-at-element
%   VERTEXAGE  vertex age
%   ELEMENTAGE element age
%
% Author: M. Moller, TU Delft, 2014.

narginchk(1,inf);
nargoutchk(nargin-1,nargin-1);

if ~all(cellfun(@ischar,varargin))
    error('Not all input arguments are strings.')
end

global lp_meshadapt

% Get primitive data from adaptation structure
nel = calllib('meshadapt', 'madapt_getnel', lp_meshadapt, useTria);
nvt = calllib('meshadapt', 'madapt_getnvt', lp_meshadapt, useTria);
ndim = calllib('meshadapt', 'madapt_getndim', lp_meshadapt, useTria);
nnve = calllib('meshadapt', 'madapt_getnnve', lp_meshadapt, useTria);


for i=1:length(varargin)
    switch lower(varargin{i})
      case 'nel'
        varargout{i} = nel;
      case 'nvt'
        varargout{i} = nvt;
      case 'ndim'
        varargout{i} = ndim;
      case 'nnve'
        varargout{i} = nnve;
      case 'coords'
        p_coords = libpointer('doublePtr', zeros(ndim,nvt));
        calllib('meshadapt', 'madapt_getvertexcoords', ...
                lp_meshadapt, p_coords, useTria);
        varargout{i} = p_coords.Value;
      case 'vertices'
        p_vertices = libpointer('int32Ptr', zeros(nnve,nel));
        calllib('meshadapt', 'madapt_getverticesatelement', ...
                lp_meshadapt, p_vertices, useTria);
        varargout{i} = p_vertices.Value;
      case 'neighbours'
        p_neighbours = libpointer('int32Ptr', zeros(nnve,nel));
        calllib('meshadapt', 'madapt_getneighboursatelement', ...
                lp_meshadapt, p_neighbours, useTria);
        varargout{i} = p_neighbours.Value;
      case 'vertexage'
        p_vertexage = libpointer('int32Ptr', zeros(1,nvt));
        calllib('meshadapt', 'madapt_getvertexage', ...
                lp_meshadapt, p_vertexage);
        varargout{i} = p_vertexage.Value;
      case 'elementage'
        p_elementage = libpointer('int32Ptr', zeros(1,nel));
        calllib('meshadapt', 'madapt_getelementage', ...
                lp_meshadapt, p_elementage);
        varargout{i} = p_elementage.Value;
      otherwise
        error(['Invalid keyword' upper(varargin{i})])
    end 
end
