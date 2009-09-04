function animate

load out.dat
load parameter.dat

nx = parameter;

nstep=size(out,1)/nx;

x=reshape(out(:,1),nx,nstep);
h=reshape(out(:,2),nx,nstep);
u=reshape(out(:,3),nx,nstep);

towait=0.0;

% Animation
pause on;
% eigentlich den hier: plot(x(:,1),h(:,1),'.-',x(:,1),u(:,1),'.-')
plot(x(:,1),h(:,1),'.-')
  %set(gca,'XLim',[xl xr]);
  %set(gca,'YLim',[-0.5 1.5]);
  drawnow;
pause
for n=1:10:nstep
  clf
% eigentlich den hier: plot(x(:,1),h(:,1),'.-',x(:,1),u(:,1),'.-')
plot(x(:,1),h(:,1),'.-')
  %set(gca,'XLim',[xl xr]);
  %set(gca,'YLim',[-0.5 1.5]);
  drawnow;
  pause(towait);
end   