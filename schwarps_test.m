% schwarps example
addpath(genpath('BBS'));
addpath(genpath('schwarps'));

er = 1e-4;
t= 1e-3;
nC = 20; % number of control centers of the BBS 
par = 1e-3; %schwarzian parameter
p = 20; % size of grid

load data_schwarps % q1 = points on first image, q2 = points on second image
% Goal : to estimate a warp between image 1 and 2
umin = min(q1(1,:))-t; umax = max(q1(1,:))+t;
vmin = min(q1(2,:))-t; vmax = max(q1(2,:))+t;
bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 2);
coloc = bbs_coloc(bbs, q1(1,:), q1(2,:));
lambdas = er*ones(nC-3, nC-3);
bending = bbs_bending(bbs, lambdas);
[xv,yv]=meshgrid(linspace(bbs.umin,bbs.umax,p),linspace(bbs.vmin,bbs.vmax,p));

% grid on first image
g1(1,:) = xv(:)';
g1(2,:) = yv(:)';

cpts = (coloc'*coloc + bending) \ (coloc'*q2(1:2,:)');
ctrlpts = cpts';

% schwarps correction
ctrlpts = optimPanalSchwarz(bbs,ctrlpts,q1(1:2,:)',q2(1:2,:)',[xv(:),yv(:)],par);
ctrlpts =ctrlpts';
 
q2_est =  bbs_eval(bbs,ctrlpts,q1(1,:)',q1(2,:)',0,0);
g2_est =  bbs_eval(bbs,ctrlpts,g1(1,:)',g1(2,:)',0,0);

error=sqrt(mean((q2_est(1,:)-q2(1,:)).^2+(q2_est(2,:)-q2(2,:)).^2));
disp(fprintf('[ETA] Warp estimation error error = %f',error));

figure(1)
title('Grid and points on image 1')
hold on;
plot(q1(1,:),q1(2,:),'ro');
mesh(reshape(g1(1,:),size(xv)),reshape(g1(2,:),size(xv)),zeros(size(xv)),'FaceColor','none');
axis equal
hold off;

figure(2)
title('Warped Grid and points on image 2')
hold on;
plot(q2(1,:),q2(2,:),'ro');
plot(q2_est(1,:),q2_est(2,:),'g*');
mesh(reshape(g2_est(1,:),size(xv)),reshape(g2_est(2,:),size(xv)),zeros(size(xv)),'FaceColor','none');
axis equal
hold off;
