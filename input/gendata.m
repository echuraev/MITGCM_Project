% This is a matlab script that generates the input data

% Questions:
% What is the number of Reynolds?
%

clear prec ieee dx x dy y dz z i j k t grid_x grid_z zp s fid eta r0 ...
    r_m d r0 rc rloc_max rloc_min tAlpha zpm zpd CWD h par_path ans;

% precision for output files
prec='real*8';
ieee='b';

%** Domain dimmension parameters
%*  begin
% Area parameters reading
par_path = './shares';
CWD = pwd;
cd(par_path);
h = parameters;
cd(CWD);
%* end

%** X resolution (m)
%*  begin
dx=zeros(h.nx,1);
x=zeros(h.nx,1);
for i=1:h.nx
	dx(i) = h.Lx/(h.nx);
end;
x(1) = 0;
for i=2:h.nx
    x(i) = x(i - 1) + dx(i);
end;
%** end

%** Y resolution (m)
%*  begin
dy=zeros(h.ny,1);
y=zeros(h.ny,1);
if (h.ny ~= 1)
    for j=1:h.ny
        dy(j) = h.Ly/(h.ny);
    end;
else
    dy(1) = h.Lx/(h.nx+1);
end;
y(1) = dy(1);
for j=2:h.ny
    y(j) = y(j - 1) + dy(j);
end;
%** end

%** R resolution (m)
%*  begin
dz=h.H*ones(h.nz,1)/h.nz;
z=-dz(1)/2:-dz(1):-h.H; % Cell centered points
%** end

%** Temperature profile
%*  begin
zpm = (h.zpr + h.zpl)/2;
zpd = zpm - h.zpr;
t = zeros(h.nx, h.nz);
for i=1:h.nx
    for k=1:h.nz
        t(i, k) = h.dr*tanh((z(k) - zpm + zpd*tanh((x(i) - h.Lg)/h.dg))/h.dh);
        if (z(k) > (h.zpl + h.dh) && z(k) < (h.zpr - h.dh))
            t(i, k) = -h.dr*tanh((x(i) - h.Lg)/h.dg);
        end;
    end;
end;
if (h.deb ~= 0)
    figure();
    grid on;
    [grid_x, grid_z] = meshgrid(x, z);
    pcolor(grid_x', grid_z', t);
    shading interp;
    colormap(jet);
    colorbar;
    title(sprintf('Temperature'));
end
fid=fopen('temperature.init','w',ieee); fwrite(fid,t,prec); fclose(fid);
%** end.

%** Salinity profile
%*  begin
s(h.nx, h.ny, h.nz)=0.0;%*rand([nx,ny,nz]);
for k=1:h.nz
    s(:,:,k) = 0.0;
end;
fid=fopen('salinity.init','w',ieee); fwrite(fid,s,prec); fclose(fid);
%** end.

%** Sloping channel
%*  begin
d = zeros(h.nx, h.ny);
for i=1:h.nx
    for j=1:h.ny
        d(i,j) = -h.H;
    end;
    if (i == h.nx || i == 1)
        d(i, j) = 0;
    end;
end;
if (h.deb ~= 0)
    figure();
    plot(x, d(:,1));
    title(sprintf('Topography'));
end;
fid=fopen('topography.init','w',ieee); fwrite(fid,d,prec); fclose(fid);
%** end

%** Surface elevation
%** begin
eta = zeros(h.nx, h.ny);
for i=1:h.nx
    for j=1:h.ny
        if (x(i) < h.Lg)
            eta(i, j) = 0;
        else
            eta(i, j) = -h.dSurf;
        end;
    end;
end;
if (h.deb ~= 0)
    figure();
    plot(x, eta(:,1));
    title(sprintf('Surface'));
end;
fid=fopen('surface.init','w',ieee); fwrite(fid,eta,prec); fclose(fid);
%** end

%** Save coordinate steps arrays
%* begin
fid=fopen('delXvar', 'w', ieee); fwrite(fid, dx, prec); fclose(fid);
fid=fopen('delYvar', 'w', ieee); fwrite(fid, dy, prec); fclose(fid);
fid=fopen('delRvar', 'w', ieee); fwrite(fid, dz, prec); fclose(fid);
%* end

%** Print variables
%* begin
disp('    Parameters, which were used in article:');
disp(['1. Density of upper layer r2: 1022, Density of lower layer r1: 1047']);
disp(['2. Depth of added water column hv: 0.158']);
disp(['3. Main section increase dH: 0.0143']);
disp(['4. Total depth H: 0.265, depth of the left pyc h2*: 0.21, depth of the right pyc h2: 0.05']);
disp(['5. Internal long wave velocity c0: 0.097']);
disp(['6. Reynolds number Rew: 2.6*10e4']);

disp('    Parameters, which are used in simulation:    ');
disp(['1. Density of upper layer r2: ' num2str(h.r2) ', Density of lower layer r1: ' num2str(h.r1)]);
disp(['2. Depth of added water column hv: ' num2str(h.h2add, 8)]);
disp(['3. Main section increase dH: ' num2str(h.H - 0.25, 8)]);
disp(['4. Total depth H: ' num2str(h.H, 8) ', depth of the left pyc h2*: ' num2str(h.h2l, 8) ...
    ', depth of the right pyc h2: ' num2str(h.h2r, 8)]);
disp(['5. Internal long wave velocity c0: ' num2str(h.c0, 8)]);
disp(['6. Reynolds number Rew: ' num2str(h.ReW, 8)]);
disp(['    Calculated density:']);
tAlpha = 9.6664e-4;
r0 = 1034.5;
rc = 1034.5;
rloc_max = r0*(-min(min(t))*tAlpha);
rloc_min = r0*(-max(max(t))*tAlpha);
disp(['Reference density: lower ' num2str(rloc_max, 8) ' upper ' num2str(rloc_min, 8)]);
disp(['Total density: lower ' num2str(rc + rloc_max, 8) ' upper ' num2str(rc + rloc_min, 8)]);
%* end

clear prec ieee dx x dy y dz z i j k t grid_x grid_z zp s fid eta r0 ...
    r_m d r0 rc rloc_max rloc_min tAlpha zpm zpd CWD h par_path ans;
