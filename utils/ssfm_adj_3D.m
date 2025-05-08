function [dmig] = ssfm_adj_3D(din,vel_p,vel_s,fpeak,timelag,nt,dt,dx,dy,dz,f1,f2,bc,src_type,src_func,save_wavefield,if_cg)

% SSA 偏移
if ~exist('if_cg','var')
  if_cg = 0;
end

[nz,nx,ny] = size(vel_p);

if if_cg==1
   din = reshape(din,nt,nx,ny);
end

dmig = zeros(nz,nx,ny);

% 频率点数
nf = 2^nextpow2(nt);
w = 2*pi/(nf*dt)*[0:(nf/2-1),-(nf/2):-1];
ind= (w==0);
w(ind)=1E-30;

%define loop over frequencies
iw1 = floor(f1*dt*nf)+1;
iw2 = floor(f2*dt*nf)+1;

if iw2 > floor(nf/2)+1
    iw2=floor(nf/2)+1;
end


%% 设置平均速度和扰动
du_p = zeros(nz,nx,ny);  
vavg_p = zeros(nz,1);
du_s = zeros(nz,nx,ny);
vavg_s = zeros(nz,1);

for iz = 1:nz
    vp = squeeze(vel_p(iz,:,:));
    vavg_p(iz) = mean(vp(:));             % 每一层平均速度 纵波
    vs = squeeze(vel_s(iz,:,:));
    vavg_s(iz) = mean(vs(:));             % 每一层平均速度 横波
    for ix = 1:nx
        for iy = 1:ny
            du_p(iz,ix,iy) = (1/vel_p(iz,ix,iy)-1/vavg_p(iz));            % 扰动速度 纵波
            du_s(iz,ix,iy) = (1/vel_s(iz,ix,iy)-1/vavg_s(iz));            % 扰动速度 横波
        end
    end
end

%%      
%---------------------------任意平面入射--------------------------------%
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

dsc = fft(dsrc,nf,1); 
% (w,x,y)
dfx = fft(din(:,:,:),nf,1);

img = zeros(nz,nx,ny);

parfor iw = iw1:iw2
    % 震源波场正传
    source_input = dsc(iw,:,:);
    source_input = reshape(source_input,nx,ny);
    [swave,~] = sspropog_op_3D(source_input,vavg_p,du_p,nx,dx,nz,dz,ny,dy,w(iw), 1,bc,'source',  save_wavefield);
    % swave(nz,nx,ny) F-X domain in frequency slice
    
    % real data
    receiver_input = reshape(dfx(iw,:,:),nx,ny);
    [rwave,~] = sspropog_op_3D(receiver_input,vavg_s,du_s,nx,dx,nz,dz,ny,dy,w(iw),-1,bc,'receiver',save_wavefield);
    % (iw,nz,nx,ny)
    %apply imaging condition
    img = img + real(rwave.*conj(swave));
    
end
dmig(:,:,:) = img./nf;        

if if_cg==1
dmig = dmig(:); 
end

return

