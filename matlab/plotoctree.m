function plotoctree(cfilename)
%
% PLOTOCTREE  Visualize an octree structure.
%
%   PLOTOCTREE(CFILENAME) reads the octree structure from the given
%     file and visualizes all vertices in the (x,y,z)-space and the cubes.

fid=fopen(cfilename,'r');

% The first entry must be a hexahedra which is used as bounding box
key=fscanf(fid,'%s\n',1);
if key=='hex',
  bbox=fscanf(fid,'%g',6);
  clf
  axis([bbox(1),bbox(4),bbox(2),bbox(5),bbox(3),bbox(6)]);
  hold on
else
  error('Missing bounding box in file');
end

% Read hexahedras and points
while 1,

  key=fscanf(fid,'%s\n',1);
  
  if length(key)==0,
    break;
  end
    
  switch lower(key)
      case {'hex'}
          hex=fscanf(fid,'%g',6);
          line(hex([1,4]),hex([2,2]),hex([3,3]));
          line(hex([1,1]),hex([2,5]),hex([3,3]));
          line(hex([1,1]),hex([2,2]),hex([3,6]));
          line(hex([4,4]),hex([2,5]),hex([3,3]));
          line(hex([4,4]),hex([2,2]),hex([3,6]));
          line(hex([1,4]),hex([5,5]),hex([3,3]));
          line(hex([1,4]),hex([5,5]),hex([6,6]));
          line(hex([1,4]),hex([2,2]),hex([6,6]));
          line(hex([4,4]),hex([2,5]),hex([6,6]));
          line(hex([1,1]),hex([2,5]),hex([6,6]));
          line(hex([1,1]),hex([5,5]),hex([3,6]));
          line(hex([4,4]),hex([5,5]),hex([3,6]));
          
      case {'node'}
          node=fscanf(fid,'%g',3);
          plot3(node(1),node(2),node(3),'.r');
  end
  
end
fclose(fid);


