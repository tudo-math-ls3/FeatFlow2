function plotquadtree(cfilename)
%
% PLOTQUADTREE  Visualize a quadtree structure.
%
%   PLOTQUADTREE(CFILENAME) reads the quadtree structure from the given
%     file and visualizes all vertices in the (x,y)- plane and the quads.

fid=fopen(cfilename,'r');

% The first entry must be a rectangle which is used as bounding box
key=fscanf(fid,'%s\n',1);
if key=='rect',
  bbox=fscanf(fid,'%g',4);
  clf
  axis([bbox(1),bbox(3),bbox(2),bbox(4)]);
  hold on
else
  error('Missing bounding box in file');
end

% Read rectangles and points
while 1,

  key=fscanf(fid,'%s\n',1);
  
  if length(key)==0,
    break;
  end

  switch lower(key)
      case{'rect'}
          rect=fscanf(fid,'%g',5);
          rectangle('Position',[rect(1) rect(2) rect(3)-rect(1) rect(4)-rect(2)]);
          text(0.5*(rect(1)+rect(3)), 0.5*(rect(2)+rect(4)), int2str(rect(5)));
                    
      case{'node'}
          node=fscanf(fid,'%g',2);
          plot(node(1),node(2),'.r');
  end
  
end
fclose(fid);


