function makepic(filename)

% Import numerical data
N=importdata(filename,',',1);

% Determine number of grid points
n=size(N.data,1);

% Import exact data
switch n,
    case 101
        E=importdata('exact-nlev2.dat');
    case 201
        E=importdata('exact-nlev3.dat');
    case 401
        E=importdata('exact-nlev4.dat');
    case 801
        E=importdata('exact-nlev5.dat');
    case 1601
        E=importdata('exact-nlev6.dat');
    case 3201
        E=importdata('exact-nlev7.dat');
    otherwise
        error('Unsupported number of grid points');
end        

% Extract coordinates
x=getvar(N,'Points:0');
[x,idx]=sort(x);
xE=E(:,1);

% Extract density
rho=getvar(N,'density'); rho=rho(idx);
rhoE=E(:,2);

% Extract pressure
p=getvar(N,'pressure'); p=p(idx);
pE=E(:,4);


% Extract x-velocity
u=getvar(N,'velocity:0'); u=u(idx);
uE=E(:,3);

% Plot exact and numerical solutions
plot(xE,rhoE,'--',xE,uE,'--',xE,pE,'--')
hold on
plot(x,rho,'.',x,u,'.',x,p,'.')
axis([0 1 -0.19 1.19]);
hold off

% Export figure to file
h=gcf;
[dir,name,ext] = fileparts(filename);

% Black and white
exportfig(h,[dir '/' name '-bw.eps'],...
    'Color','bw',...
    'FontMode','fixed','FontSize','16',...
    'LineMode','fixed','LineWidth','0.5');

% RGB color
exportfig(h,[dir '/' name '-rgb.eps'],...
    'Color','rgb',...
    'FontMode','fixed','FontSize','16',...
    'LineMode','fixed','LineWidth','0.5');

%CMYK color
exportfig(h,[dir '/' name '-cmyk.eps'],...
    'Color','cmyk',...
    'FontMode','fixed','FontSize','16',...
    'LineMode','fixed','LineWidth','0.5');

function v = getvar(S,varname)

for i=1:size(S.data,2)
    if ~isempty(findstr(S.textdata{i},varname));
        v=S.data(:,i);
        return
    end
end
