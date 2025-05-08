% 3D receiver function Least Squares Migration
% wirtten by Pengfei Zuo
% Zhejiang university
% Synthetic data of undulated moho
% Dec 28, 2024
clc;
clear;
close all;
addpath(genpath('./utils/'))

%% step 0: load receiver functions data
mig_undulated = struct;
disp('Loading the data...')

load('./matfiles/RF_undulated_syn.mat')  % only including receiver function data
disp('done!')
[nbaz,ninc] = size(RF_undulated_syn);

inc_angles = 15:5:35;
baz_angles = 30:30:360;
% station's location infomation
filename = './doc/STATIONS_LOCATION';
fileID = fopen(filename);
refile = textscan(fileID,'%s %s %f %f %f %f');
ry = [refile{1,3}]/1000;
rx = [refile{1,4}]/1000;
rx = unique(rx);  % in-line (km)
ry = unique(ry);  % cross-line (km)

snx = 41;
sny = 21;
clear refile

%% step 1: create regular mesh grid
dx = 5;
dy = 5;
dz = 1;
zmax = 100;
x = 0 : dx : 300;
y = 0 : dy : 300;
z = 0 : dz : zmax;

nx = length(x);
ny = length(y);
nz = length(z);

dt = TIME(2) - TIME(1);
nt = length(TIME);
%% step 2: velocity model (3D model)
npt = 81;
xl = [0,300];
yl = [0,200];

dxi = (xl(2)-xl(1))/(npt-1);
dyi = (yl(2)-yl(1))/(npt-1);
xi = xl(1) + (0:npt - 1)*dxi;
yi = yl(1) + (0:npt - 1)*dyi;

zdepth = 2*peaks(npt)+40;
zdepth = zdepth';
[xv,yv] = meshgrid(xi,yi);

[xvq,yvq] = meshgrid(x,y);
zdepth_inq = interp2(xv,yv,zdepth,xvq,yvq);
zdepth_inq(isnan(zdepth_inq)) = 40;
% figure
% surf(xvq,yvq,zdepth_inq)
% colormap(gray);axis equal
% shading interp
% set(gca,'ZDir','reverse','FontSize',26);
% xlabel('Distance (km)');
% ylabel('Distance (km)');
% zlabel('Depth (km)');
% ylim([0 200])


Vp = zeros(nz,nx);
Vs = zeros(nz,nx);

vel_p = zeros(nz,nx,ny);
vel_s = vel_p;
for i = 1:nz
    mask = ones(ny,nx)*(i-1)*dz <= zdepth_inq;
    vel_p(i,:,:) = mask' * 5800 + ~mask' * 8289; 
    vel_s(i,:,:) = mask' * 3460 + ~mask' * 4517; 
end

% smoothing
N = 5;
for i = 1:ny
    [vel_p(:,:,i),~] = moving_avg(vel_p(:,:,i),N,'constant',1);
    [vel_p(:,:,i),~] = moving_avg(vel_p(:,:,i),N,'constant',2);

    [vel_s(:,:,i),~] = moving_avg(vel_s(:,:,i),N,'constant',1);
    [vel_s(:,:,i),~] = moving_avg(vel_s(:,:,i),N,'constant',2);
end
for i = 1:nz
    vtp = reshape(vel_p(i,:,:),nx,ny);
    [vel_p(i,:,:),~] = moving_avg(vtp,N,'constant',2);

    vstp = reshape(vel_s(i,:,:),nx,ny);
    [vel_s(i,:,:),~] = moving_avg(vstp,N,'constant',2);
end
vel_p = vel_p/1000;
vel_s = vel_s/1000;


%% -------------------source parameters------------------------------------
fpeak = 1.2;        % ricker
flow = 0.01;
fhigh = 2.5;
bc = 1; 
%% step 3
for ishot = 9
    [iinc,ibaz] = ind2sub([ninc,nbaz],ishot);

    disp(['Baz: ',num2str(baz_angles(ibaz))])
    disp(['Inc: ',num2str(inc_angles(iinc))])
    take_off = inc_angles(iinc);
    back_az = baz_angles(ibaz);

    rf0 = RF_undulated_syn(ibaz,iinc).rf;
    % normaliztion
    itr = rf0./max(rf0(:));

    % calculate time shift on the bottom of the model
    % to simulate plane wave 
    v_mean = mean(vel_p(end,:,:),'all');
    timelag = calculate_timeshift(take_off,back_az,x,y,v_mean);
    % figure
    % imagesc(timelag)
    % extract source time function from RFs
    % [pcaData,projectionVectors,eigVal] = myPCA(itr,1);
    srctmp=mean(itr,2);
    [win] = filt_win(srctmp,TIME,-2,1,-1);
    src=srctmp.*win;
    src_func = src./max(src);

    % taper RF to remove later conversions
    [win] = filt_win(srctmp,TIME,1,9,3);
    win = win*ones(1,size(itr,2));
    itr=itr.*win;

    %% forward modeling
    save_wavefield = 1;
    disp(['shot : ',num2str(ishot),' forward modeling ...'])
    img = zeros(nz,nx,ny);
    tic
    [d,mod_source,~] = ssfm_fd_3D(img,vel_p,vel_s,fpeak,timelag,nt,dt,dx,dz,dy,flow,fhigh,bc,'P',src_func,save_wavefield);
    disp('done!')
    toc;

    Zstrace = reshape(mod_source(1,:,:,:),nt,nx,ny);
    % cross correlate P wave with RF
    tdelay = zeros(nx,ny);
    dshift = zeros(nt,nx,ny);

    [tq,xq,yq]=meshgrid(x,TIME,y);
    [t1,x1,y1] = meshgrid(rx,TIME,ry);
    vv=reshape(itr,nt,snx,sny);
    rftmp = interp3(t1,x1,y1,vv,tq,xq,yq,'nearest');

    it = (0:nt-1)*dt;
    for i = 1:nx
        for j = 1:ny
            temp1 = src;
            temp2 = Zstrace(:,i,j);
            xc=xcorr(temp2,temp1);
    
            tax = [-(nt-1):(nt-1)]*dt;
            [~,ind]=max(xc);
    
            tdelay(i,j) = tax(ind);

            d(:,i,j) = rftmp(:,i,j);

            dshift(:,i,j) = fftShift(d(:,i,j),it,tdelay(i,j));
        end
    end
    dshift(isnan(dshift))=0;

    %% do migration
    save_wavefield = 0;
    bc=1;
    disp('==========>')
    disp('do migration')
    tic
    [mig] = ssfm_adj_3D(dshift,vel_p,vel_s,fpeak,timelag,nt,dt,dx,dy,dz,flow,fhigh,bc,'P',src_func,save_wavefield);
    toc
    % [dp,mod_source,mod_receiver] = ssfm_fd_3D(mig,vel_p,vel_s,fpeak,timelag,nt,dt,dx,dz,dy,flow,fhigh,bc,'P',src_func,save_wavefield);
    % disp('done!')

    mig_undulated(ishot).mig= mig;

    % S=repmat(squeeze(any(dshift)),1,1,size(dshift,1));
    % S=permute(S,[3,1,2]);
    % S=S(:);
    % itermax = 10;
    % if_cg = 1;
    % L = @(m) S.*ssfm_fd_3D(m,vel_p,vel_s,fpeak,timelag,nt,dt,dx,dy,dz,flow,fhigh,bc,'P',src_func,save_wavefield,if_cg);
    % Lt= @(d) ssfm_adj_3D(S.*d,vel_p,vel_s,fpeak,timelag,nt,dt,dx,dz,dy,flow,fhigh,bc,'P',src_func,save_wavefield,if_cg);
    % m1=randn(nz,nx,ny);
    % d2=randn(nt,nx,ny);
    % m1=m1(:);
    % d2=d2(:);
    % d1=L(m1);
    % m2= Lt(d2);
    % dot1=sum(sum(sum(d1.*d2)))
    % dot2=sum(sum(sum(m1.*m2)))

    %% LS with model constrain
    disp(['==========>',num2str(ishot)])
    disp('Preconditioned least-squares migration begins')
    din = dshift(:);
    S = repmat(any(dshift),size(dshift,1),1);
    S = S(:);
    mu = 0.1;
    if_cg = 1;
    P = @(u) precon_x_3d(u,nx,ny,nz,dx,dy,dz);
    Pt= @(d) preconT_x_3d(d,nx,ny,nz,dx,dy,dz);

    L = @(m) S.*ssfm_fd_3D(m,vel_p,vel_s,fpeak,timelag,...
        nt,dt,dx,dz,dy,flow,fhigh,bc,'P',src_func,save_wavefield,if_cg);
    
    Lt = @(d) ssfm_adj_3D(S.*d,vel_p,vel_s,fpeak,timelag,...
        nt,dt,dx,dy,dz,flow,fhigh,bc,'P',src_func,save_wavefield,if_cg);    
    LP = @(u) L(P(u));
    PtLt = @(d) Pt(Lt(d));

    %-----------------------------
    % % dot product test
    % m1=randn(nz,nx,ny);
    % d2=randn(nt,nx,ny);
    % m1=m1(:);
    % d2=d2(:);
    % 
    % d1=LP(m1);
    % m2= PtLt(d2);
    % 
    % dot1=sum(sum(sum(d1.*d2)));
    % dot2=sum(sum(sum(m1.*m2)));
    % 
    % fprintf('dot1: %3f \n',dot1);
    % fprintf('dot2: %3f \n',dot2);
    % if abs(dot1 - dot2) < 10^-10
    %     disp('pass the dot product  test! ')
    % else
    %     disp('try again...')
    % end
    %--------------------------------

    A = @(u) PtLt(LP(u)) + mu * u;
    b = PtLt(din);

    tic;
    itermax = 5;
    disp('PCG begins...')
    [utmp1,flag,relres,iter,resvec] = pcg(A,b,[],itermax);
    disp(['flag:',num2str(flag),', min value:',num2str(min(resvec))])
    toc;

    % migration image
    
    mtmp1=P(utmp1(:));
    mtmp1=reshape(mtmp1,nz,nx,ny);
    % save data
    mig_undulated(ishot).lsmp = mtmp1;

    % predict data
    % Hd1 = L(mtmp1(:));
    % pre_dp1 = reshape(Hd1,nt,nx,ny);

    % image_step(ishot).prelsmp = pre_dp1;

    nnx = find(y==100);    % y axis  km
    nny = find(x==150);    % x axis  km
    figure
    set(gcf,'Position',[1,200,1950,700],'Color','white')
    subplot(2,5,[1 2 3])
    imagesc(x,z,mig(:,:,nnx)./max(mig(:)),'CDataMapping','scaled','Interpolation','bilinear');
    title('Migration image')
    clim([-1 1]);colormap('jet')
    axis equal
    ylim([0 100])
    xlim([0 300])
    clim([-1 1])
    colorbar
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',16)
    xlabel('X (km)')
    ylabel('Depth (km)')
    hold on
    scatter(rx,zeros(1,snx),30,'v','filled',...
                  'MarkerFaceColor',[1 0 0]);
    plot(x,zdepth_inq(nnx,:),'--','LineWidth',1.5,'Color','black')
    hold off
    text(-15,-12,'(a)','FontSize',20)
    
    
    subplot(2,5,[6 7 8])
    imagesc(x,z,mtmp1(:,:,nnx)./max(mtmp1(:)),'CDataMapping','scaled','Interpolation','bilinear');
    title('LSM image')
    axis equal
    ylim([0 100])
    xlim([0 300])
    clim([-1 1]);colormap('jet')
    colorbar
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',16)
    xlabel('X (km)')
    ylabel('Depth (km)')
    hold on
    scatter(rx,zeros(1,snx),30,'v','filled',...
                  'MarkerFaceColor',[1 0 0]);
    plot(x,zdepth_inq(nnx,:),'--','LineWidth',1.5,'Color','black')
    hold off
    text(-15,-12,'(c)','FontSize',20)
    
    subplot(2,5,[4 5])
    imagesc(y,z,squeeze(mig(:,nny,:))./max(mig(:)),'CDataMapping','scaled','Interpolation','bilinear');
    title('Migration image')
    clim([-1 1]);colormap('jet')
    axis equal
    ylim([0 100])
    xlim([0 200])
    % xlim([0 300])
    clim([-1 1])
    colorbar
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',16)
    xlabel('Y (km)')
    ylabel('Depth (km)')
    hold on
    scatter(ry,zeros(1,sny),30,'v','filled',...
                  'MarkerFaceColor',[1 0 0]);
    plot(y,zdepth_inq(:,nny),'--','LineWidth',1.5,'Color','black')
    hold off
    text(-15,-12,'(b)','FontSize',20)
    
    
    subplot(2,5,[9 10])
    imagesc(y,z,squeeze(mtmp1(:,nny,:))./max(mtmp1(:)),'CDataMapping','scaled','Interpolation','bilinear');
    title('LSM image')
    axis equal;colormap('jet')
    ylim([0 100])
    % xlim([0 300])
    clim([-1 1])
    colorbar
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',16)
    xlabel('Y (km)')
    ylabel('Depth (km)')
    hold on
    scatter(ry,zeros(1,sny),30,'v','filled',...
                  'MarkerFaceColor',[1 0 0]);
    plot(y,zdepth_inq(:,nny),'--','LineWidth',1.5,'Color','black')
    hold off
    xlim([0 200])
    text(-15,-12,'(d)','FontSize',20)

    % figdir = './figures/undulated_model/';
    % figname = [figdir,'img_baz_',num2str(back_az),'_inc_',num2str(take_off),'.png'];
    % export_fig(figname)
    close all
end

%% plot stacked
mig_stacked = zeros(nz,nx,ny);
miglsp = zeros(nz,nx,ny);
for i = 1:length(mig_undulated)
    mig_stacked = mig_stacked + mig_undulated(i).mig;
    miglsp = miglsp + mig_undulated(i).lsmp;
end
mig_stacked = mig_stacked./length(mig_undulated);
miglsp = miglsp./length(mig_undulated);

nnx = find(y==100);    % y axis  km
nny = find(x==150);    % x axis  km
figure
set(gcf,'Position',[1,200,1950,700],'Color','white')
subplot(2,5,[1 2 3])
imagesc(x,z,mig_stacked(:,:,nnx)./max(mig_stacked(:)),'CDataMapping','scaled','Interpolation','bilinear');
title('Migration image (stacked)')
clim([-1 1]);
axis equal
colormap('jet')
ylim([0 100])
xlim([0 300])
clim([-1 1])
colorbar
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'FontSize',16)
xlabel('X (km)')
ylabel('Depth (km)')
hold on
scatter(rx,zeros(1,snx),30,'v','filled',...
              'MarkerFaceColor',[1 0 0]);
plot(x,zdepth_inq(nnx,:),'--','LineWidth',1.5,'Color','black')
hold off
text(-15,-12,'(a)','FontSize',20)

subplot(2,5,[6 7 8])
imagesc(x,z,miglsp(:,:,nnx)./max(miglsp(:)),'CDataMapping','scaled','Interpolation','bilinear');
title('LSM image (stacked)')
axis equal
colormap('jet')
ylim([0 100])
xlim([0 300])
clim([-1 1])
colorbar
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'FontSize',16)
xlabel('X (km)')
ylabel('Depth (km)')
hold on
scatter(rx,zeros(1,snx),30,'v','filled',...
              'MarkerFaceColor',[1 0 0]);
plot(x,zdepth_inq(nnx,:),'--','LineWidth',1.5,'Color','black')
hold off
text(-15,-12,'(c)','FontSize',20)

subplot(2,5,[4 5])
imagesc(y,z,squeeze(mig_stacked(:,nny,:))./max(mig_stacked(:)),'CDataMapping','scaled','Interpolation','bilinear');
title('Migration image (stacked)')
clim([-1 1]);
axis equal
colormap('jet')
ylim([0 100])
xlim([0 200])
% xlim([0 300])
clim([-1 1])
colorbar
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'FontSize',16)
xlabel('Y (km)')
ylabel('Depth (km)')
hold on
scatter(ry,zeros(1,sny),30,'v','filled',...
              'MarkerFaceColor',[1 0 0]);
plot(y,zdepth_inq(:,nny),'--','LineWidth',1.5,'Color','black')
hold off
text(-15,-12,'(b)','FontSize',20)

subplot(2,5,[9 10])
imagesc(y,z,squeeze(miglsp(:,nny,:))./max(miglsp(:)),'CDataMapping','scaled','Interpolation','bilinear');
title('LSM image (stacked)')
axis equal
colormap('jet')
ylim([0 100])
% xlim([0 300])
clim([-1 1])
colorbar
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'FontSize',16)
xlabel('Y (km)')
ylabel('Depth (km)')
hold on
scatter(ry,zeros(1,sny),30,'v','filled',...
              'MarkerFaceColor',[1 0 0]);
plot(y,zdepth_inq(:,nny),'--','LineWidth',1.5,'Color','black')
hold off
xlim([0 200])
text(-15,-12,'(d)','FontSize',20)

% figdir = './figures/undulated_model/';
% figname = [figdir,'img_undulated_stacked.png'];
% export_fig(figname)

%% plot RF
figure
set(gcf,'Position',[200 300 1250 550],'Color','white')
subplot(121)
wigb(squeeze(dshift(:,nny,:)),1,y,TIME)
xlabel('Distance (km)')
ylabel('Time (sec)')
title('Shifted RFs (Ps phase)')
set(gca,'FontSize',18)
xlim([0 200])
subplot(122)
wigb(squeeze(Zstrace(:,nny,:)),1,y,TIME)
xlabel('Distance (km)')
ylabel('Time (sec)')
title('Records')
set(gca,'FontSize',18)
xlim([0 200])

% 
% figname = 'undulate_miglsp_3d.png';
% export_fig(figname)

