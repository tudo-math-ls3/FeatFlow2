function verify(MC, ML, A, T, Tbdr)

deg=zeros(size(MC,1),1);

[I,J]=find(MC);
for idx=1:nnz(I),
    deg(I(idx)) = deg(I(idx))+1;
end

bar(deg);

% Statistical output
display(['NEQ   = ' num2str(size(A,1))]);
display(['NCOLS = ' num2str(size(A,2))]);
display(['NA    = ' num2str(nnz(A))]);
display(['DEG   = ' num2str(max(deg))]);

bMC=0; minMC=0;

% Check sign of mass matrix entries
[I,J]=find(MC<0);

if ~isempty(I),
    for idx=1:nnz(I),
        minMC=min([minMC, MC(I(idx),J(idx))]);
    end
    disp(['Consistent mass matrix has negative entries minval='...
          num2str(minMC,4)]);
    bMC=1;
end
clear I J;


bML=0; minML=0;

% Check sign of lumped mass matrix
[I,J]=find(ML<0);

if ~isempty(I),
    for idx=1:nnz(I),
        minML=min([minML, ML(I(idx),J(idx))]);
    end
    disp(['Lumped mass matrix has negative entries minval='...
          num2str(minML,4)]);
    bML=1;
end
clear I J;


bT=0; minT=0;

% Check sign of transport operator
[I,J]=find(T);

for idx=1:nnz(I),
    i=I(idx); j=J(idx);
    % Skip diagonal entries
    if (i==j), continue, end
    % Check sign of off-diagonal entries
    if (T(i,j)<0),
        minT=min([minT, T(i,j)]);
        bT=1;
    end
end
clear I J;
if (bT~=0),
    disp(['Transport operator has negative off-diagonal entries minval='...
          num2str(minT,4)]);
end

% Status report
if (bMC==0)
    disp('Consistent mass matrix:   OK!');
end

% Status report
if (bML==0)
    disp('Lumped mass matrix:       OK!');
end

% Status report
if (bT==0)
    disp('Transport operator:       OK!');
end


