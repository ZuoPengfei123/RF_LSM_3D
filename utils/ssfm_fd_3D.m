function [mod,mod_source,mod_receiver] = ssfm_fd_3D(img,vel_p,vel_s,fpeak,timelag,nt,dt,dx,dz,dy,f1,f2,bc,src_type,src_func,save_wavefield,if_cg)

if ~exist('if_cg','var')
  if_cg = 0;
end
[nz,nx,ny] = size(vel_p);

if if_cg==1
   img = reshape(img,nz,nx,ny);
end

% receiver field
mod1 = zeros(nt,nx,ny);
nf  = 2^nextpow2(nt);
nsx = 2^nextpow2(nx);
nsy = 2^nextpow2(ny);

w = 2*pi/(nf*dt)*[0:(nf/2-1),-(nf/2):-1]; 

ind = (w==0);
w(ind)=1E-30;

iw1 = floor(f1*dt*nf)+1;
iw2 = floor(f2*dt*nf)+1;
% if frequency exceeding the nyquest freq
if iw2 > floor(nf/2)+1
    iw2 = floor(nf/2)+1;
end

%% 设置平均速度和扰动
du_p = zeros(nz,nx,ny);
vavg_p = zeros(nz,1);           
du_s = zeros(nz,nx,ny);
vavg_s = zeros(nz,1);

for iz = 1:nz
    vp = squeeze(vel_p(iz,:,:));
    vavg_p(iz) = mean(vp(:)); 
    vs = squeeze(vel_s(iz,:,:));
    vavg_s(iz) = mean(vs(:));
    for ix = 1:nx
        for iy = 1:ny
            du_p(iz,ix,iy) = (1/vel_p(iz,ix,iy) - 1/vavg_p(iz));
            du_s(iz,ix,iy) = (1/vel_s(iz,ix,iy) - 1/vavg_s(iz));

        end
    end
end

%% Source
dsrc = zeros(nt,nx,ny);
t = (0:nt-1)*dt;
for indx = 1:nx
     for indy = 1:ny         
         if strcmp(src_type,'ricker')
             par = pi*fpeak*(t-5);
             src_func = 10*exp(-par.*par).*(1-2*par.*par);  % ricker wavelet
             dsrc(:,indx,indy) = fftShift(src_func',t,timelag(indy,indx));

         else
             dsrc(:,indx,indy) = fftShift(src_func,t,timelag(indy,indx));

         end
     end
 end
% fft
% figure;
% subplot(121)
% imagesc(1:61,t,squeeze(dsrc(:,:,11)))
% subplot(122)
% imagesc(timelag)

dsc = fft(dsrc,nf,1);

% output data in frequence domian
outf = zeros(nf,nsx,nsy);
outf_source = zeros(nz,nf,nsx,nsy);
outf_receiver = zeros(nz,nf,nsx,nsy);


simg = img(:,:,:);      % (nz,nx,ny)

parfor iw = iw1:iw2
    source_input = reshape(dsc(iw,:,:),nx,ny);
    [swave,swave_full] = sspropog_op_3D(source_input,vavg_p,du_p,nx,dx,nz,dz,ny,dy,w(iw),1,bc,'source',save_wavefield);
    % swave(nz,nx,ny) F-X domain for each frequency
    % swave_full(nz,nsx,nsy) F-K domain for each frequency
    
    receiver_input = (simg.*swave);
    [rwave_s,rwave_full_s] = sspropog_op_3D(receiver_input,vavg_s,du_s,nx,dx,nz,dz,ny,dy,w(iw),1,bc,'receiver',save_wavefield);
    % rwave_s(nsx,nsy) : f-k domain of each iw  at depth zi=1
    % rwave_full_s(nz,nsx,nsy): f-k domain of each iw
    outf(iw,:,:) = rwave_s;
    
    if save_wavefield
        outf_source(:,iw,:,:) = swave_full;
        outf_receiver(:,iw,:,:) = rwave_full_s;                 
        % outf_receiver(:,iw,:)=rwave_full_s+rwave_full_p; % only used
        % when simulating P phase in RF
    end
    
end

r = real(ifft(ifft2(outf,nf,nsx),nsy,3));
%--------------transform outf to tx domain-------------

mod1(:,:,:) = r(1:nt,1:nx,1:ny);

if save_wavefield
    mod_source = zeros(nz,nt,nx,ny);
    mod_receiver = zeros(nz,nt,nx,ny); 

    parfor iz=1
        tmp_source = reshape(outf_source(iz,:,:,:),nf,nsx,nsy);
        tmp_receiver = reshape(outf_receiver(iz,:,:,:),nf,nsx,nsy);

        s_full = real(ifft(ifft2(tmp_source,nf,nsx),nsy,3));
        r_full = real(ifft(ifft2(tmp_receiver,nf,nsx),nsy,3));
        mod_source(iz,:,:,:) = s_full(1:nt,1:nx,1:ny);
        mod_receiver(iz,:,:,:) = r_full(1:nt,1:nx,1:ny);
    end
else
    mod_source=[];
    mod_receiver=[];
end

mod = mod1;   

if if_cg==1
    mod=mod1(:);
end
end
