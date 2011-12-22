function [umin umax] = plotsol(filename,lbound,ubound)

sol=load(filename);

whos

x=sol(:,1);
y=sol(:,2);
u=sol(:,3);

idx=find(u>ubound);
xu=x(idx);
yu=y(idx);
uu=u(idx);

idx=find(u<lbound);
xl=x(idx);
yl=y(idx);
ul=u(idx);

plot3(x,y,u,'ok');
hold on

if ~isempty(xl),
    plot3(xl,yl,ul,'*b');
end

if ~isempty(xu),
    plot3(xu,yu,uu,'*r');
end
hold off

umin=min(u);
umax=max(u);
