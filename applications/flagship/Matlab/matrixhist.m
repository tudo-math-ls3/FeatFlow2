function matrixhist(M)

if size(M,1) > size(M,2),
    [I J]=find(M);
else
    [J I]=find(M);
end

dof=zeros(length(M),1);

for idx=1:length(I),
    dof(I(idx)) = dof(I(idx))+1;
end

bar(dof)

dof