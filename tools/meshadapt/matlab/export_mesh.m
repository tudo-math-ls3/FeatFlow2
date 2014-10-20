% function for exporting triangulation to meshadapt compatible format
function boundary_set = export_mesh(coords,vertices,name)
% basic assumptions:
% 1) linear triangles
% 2) elements are rows of vertices
% 3) points are rows of coords
% 4) 2D only

% check for correct format input
[nI,nJ] = size(coords);
if nI==2
    coords = coords';
    n_points = nJ;
elseif nJ==2
    n_points = nI;
else
    error('Coordinates are not in 2D format');
end
[nI,nJ] = size(vertices);
if nI==3
    vertices = vertices';
    n_elems = nJ;
elseif nJ==3
    n_elems = nI;
else
    error('Triangles are not in correct format');
end

% check for CCW triangulation
X = coords(:,1);
Y = coords(:,2);
I = find(ispolycw(mat2cell(X(vertices(:,[1:3 1])),ones(n_elems,1)),mat2cell(Y(vertices(:,[1:3 1])),ones(n_elems,1))));
vertices(I,:) = fliplr(vertices(I,:));
    
TR = triangulation(vertices, coords);

% first create the information for the boundaries
% a) find entire free boundary of triangulation
free = TR.freeBoundary;
% b) create ordered & closed sets of line-elements
[NBCT,boundary_set] = order_boundary(free,coords);
% c) create per boundary the parameter lists
param = cell(NBCT,1);
for i=1:NBCT
    NP = length(boundary_set{i})-1;
    param{i} = zeros(2*NP,2);
    for j=1:NP
        param{i}(1+(j-1)*2,1) = X(boundary_set{i}(j));
        param{i}(1+(j-1)*2,2) = Y(boundary_set{i}(j));
        
        param{i}(2+(j-1)*2,1) = X(boundary_set{i}(j+1))-X(boundary_set{i}(j));
        param{i}(2+(j-1)*2,2) = Y(boundary_set{i}(j+1))-Y(boundary_set{i}(j));
    end
%     param{i}(1:2:end,1) = X(boundary_set{i}(1:end-1));
%     param{i}(1:2:end,2) = Y(boundary_set{i}(1:end-1));
%     param{i}(2:2:end,1) = diff(X(boundary_set{i}));
%     param{i}(2:2:end,2) = diff(Y(boundary_set{i}));
end
% d) on which boundary lies which node
KNPR = zeros(size(X));
for i=1:NBCT
    KNPR(boundary_set{i}) = i;
end
% e) min and max indices of boundaries
KMM = [cellfun(@min,boundary_set);cellfun(@max,boundary_set)];

% now export to the PRM file, i.e. boundary information
fid = fopen([name '.prm'],'w');

fprintf(fid,'NBCT\n','char'); %#ok<CTPCT>
fprintf(fid,'%d\n',NBCT);
for i=1:NBCT
    fprintf(fid,'IBCT\n','char'); %#ok<CTPCT>
    fprintf(fid,'%d\n',i);
    fprintf(fid,'NCOMP\n','char'); %#ok<CTPCT>
    fprintf(fid,'%d\n',length(boundary_set{i})-1);
    fprintf(fid,'ITYP NSPLINE NPAR\n','char'); %#ok<CTPCT>
    fprintf(fid,'%d %d %d\n',repmat([1 1 2],length(boundary_set{i})-1,1)');
end
fprintf(fid,'PARAMETERS\n','char'); %#ok<CTPCT>
for i=1:NBCT
    fprintf(fid,'%3.16f %3.16f\n',param{i}');
end

fclose(fid);

% convert boundary points to parametrization
for i=1:NBCT
    X(boundary_set{i}(1:end-1)) = 0:(length(boundary_set{i})-2);
    Y(boundary_set{i}(1:end-1)) = 0;
end

% now export to the TRI file, i.e. domain information
fid = fopen([name '.tri'],'w');

fprintf(fid,'Gibberish on the first line\n','char'); %#ok<CTPCT>
fprintf(fid,'Gibberish on the second line\n','char'); %#ok<CTPCT>
fprintf(fid,'%d %d 0 3 %d    NEL NVT NMT NVE NBCT\n',n_elems,n_points,NBCT);
fprintf(fid,'DCORVG\n','char'); %#ok<CTPCT>
fprintf(fid,'%3.16f %3.16f\n',[X(:) Y(:)]');
fprintf(fid,'KVERT\n','char'); %#ok<CTPCT>
fprintf(fid,'%d %d %d\n',vertices');
fprintf(fid,'KNPR\n','char'); %#ok<CTPCT>
fprintf(fid,'%d\n',KNPR);
fprintf(fid,'KMM\n','char'); %#ok<CTPCT>
fprintf(fid,'%d %d\n',KMM);

fclose(fid);

return

end

% function to create ordered & closed sets of line-elements
function [NBCT,set_out] = order_boundary(boundary,coords)

if isempty(boundary)
    NBCT = 0;
    set_out = cell(0,1);
    return
end

% start with choosing the first element as a start:
while ~isempty(boundary)
    current_set = boundary(1,:);
    boundary(1,:) = [];
    start_point = current_set(1);
    current_point = current_set(2);
    
    while start_point~=current_point
        
        [I,~] = find(boundary==current_point);
        I = I(1);
        current_set = [current_set setdiff(boundary(I,:),current_point)]; %#ok<AGROW>
        current_point = current_set(end);
        boundary(I,:) = [];
        
    end
    
    if ~exist('set','var')
        set{1} = current_set;
    else
        set{end+1} = current_set; %#ok<AGROW>
    end
    
    
end

NBCT = length(set);

% first find the enclosing boundary
K = convhull(coords);
I = find(cellfun(@(b) any(ismember(b,K)),set));
CCW = ~cellfun(@(b) ispolycw(coords(b,1),coords(b,2)),set);

set_out = cell(size(set));
counter = 0;
if any(~CCW(I))
    for i=1:length(I)
        set{I(i)} = fliplr(set{I(i)}); %#ok<AGROW>
    end
end
for i=1:length(I)
    counter = counter+1;
    set_out{counter} = set{I(counter)};
end
set(I) = [];
CCW(I) = [];

if ~isempty(set)
    % assume all other boundaries are singular holes, so must have clokwise
    % orientation
    if any(CCW)
        for i=1:length(set)
            set{i} = fliplr(set{i});
        end
    end
    set_out{(counter+1):end} = set{:};
end

end
