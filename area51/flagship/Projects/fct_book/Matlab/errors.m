function [rho_err1, p_err1, u_err1] = errors(filename)

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

% Compute lumped mass matrix
dx=(x(end)-x(1))/(n-1);
M=dx*ones(n,1); M(1)=dx/2; M(end)=dx/2;

% Compute L1-errors
rho_err1 = sum(M.*abs(rho-rhoE));
p_err1   = sum(M.*abs(p-pE));
u_err1   = sum(M.*abs(u-uE));


function v = getvar(S,varname)

for i=1:size(S.data,2)
    if ~isempty(findstr(S.textdata{i},varname));
        v=S.data(:,i);
        return
    end
end
