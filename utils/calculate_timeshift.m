function [timelag] = calculate_timeshift(take_off,back_az,x_coord,y_coord,vel)



theta = deg2rad(take_off);

back_az = back_az + 180;
phi = deg2rad(back_az);

nx = length(x_coord);
ny = length(y_coord);

timelag = zeros(nx,ny);

coord = zeros(nx*ny,2);

n1 = sin(theta)*sin(phi);
n2 = sin(theta)*cos(phi);
n3 = cos(theta);

count = 1;
for j = 1:ny
    for i = 1:nx
        coord(count,1) = x_coord(i);
        coord(count,2) = y_coord(j);
        count = count + 1;
    end
end

for k = 1:size(coord,1)
    x = coord(k,1);
    y = coord(k,2);
    sta = [x,y,0];
    kv = [n1,n2,n3];

    dist = (sta*kv')/dot([n1,n2,n3],[n1,n2,n3]);

    dtime = dist / vel;

    jj = mod(k-1,ny)+1;
    ii = fix((k-jj)/ny) + 1;

    timelag(ii,jj) = dtime;
end
if min(timelag(:))<0
    timelag = timelag + abs(min(timelag(:)));
end