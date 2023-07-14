function plotMyu(u,X,L,u_i,time_step,dt)

N = size(u,1);
h = L/N;
xgrid=zeros(N,N);
ygrid=zeros(N,N);
for j=0:(N-1)
  xgrid(j+1,:)=j*h;
  ygrid(:,j+1)=j*h;
end
pcolor(10^4*xgrid,10^4*ygrid,10^8*u(:,:,u_i))
hold on; plot(10^4*X(:,1),10^4*X(:,2),'ro')

Xplt = mod(X,L);
plot(10^4*Xplt(:,1),10^4*Xplt(:,2),'ko','MarkerSize',7,'MarkerFaceColor','k')
%             plot(10^4*Xplt(1,1),10^4*Xplt(1,2),'ro','MarkerFaceColor','r','MarkerSize',5)

axis([0,10^4*L,0,10^4*L])
hcb=colorbar;
H = ylabel(hcb, 'uM');
H.Rotation = 270;
H.Position = [4.0813   20.0000         0];
xlabel('um'); ylabel('um');
title(['t is ' num2str(dt*time_step/60,'%.1f') ' min']);
axis equal
axis manual
set(gca, 'fontsize', 18);
drawnow
