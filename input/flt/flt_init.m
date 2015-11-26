% Script allows to generate an input binary file of pkg/flt for Lagrangian
% trajectories calculating.
clear prec ieee outfilename depthfilename Lx Ly Lz Lg nx ny nz dx ...
    dy dy dz dz depth XC YC Np np Zr tstart kpart kfloat ...
    iup itop tend output fid i j k coordinates amplitude delta ind kpart_surf ...
    llength lp xcs zcs ans hpl hpr indx indzh indzl lb rb h CWD par_path ...
    Z_dim dot_p l;

%** Initial variables
%* begin
% General parameters
prec='real*8';
ieee='b';
outfilename = 'flts.init';
% Area parameters reading
par_path = '../shares';
CWD = pwd;
cd(par_path);
h = parameters;
cd(CWD);
dx = h.Lx/h.nx;
dy = h.Ly/h.ny;
dz = h.H/h.nz;
% Parameters of target arrays
xcs = 10;
zcs = 3;
lb = round(xcs/2); rb = h.nx - round(xcs/2);
Np = 0; lp = 0;
%* end

%** Center of horizontal coordinates calculating.
%* begin
XC = dx*ones(h.nx, 1);
XC = cumsum(XC, 1) - dx/2;
if (h.ny ~= 1)
    YC = dy*ones(h.ny, 1);
    YC = cumsum(YC, 1) - dy/2;
else
    YC = dx/2;
end;
Zr = -dz*ones(h.nz, 1); % There are exact number of vertical coordinates
Zr = cumsum(Zr, 1) + dz/2;
Z_dim = (h.H + Zr)/h.h2r;
%* end

%** Undersurface points calculating
%* begin
surf_p = [struct('z', Z_dim(1) - 0.1, 'ind', 0) ...
          struct('z', Z_dim(1) - 0.2, 'ind', 0) ...
          struct('z', Z_dim(1) - 0.4, 'ind', 0) ...
          struct('z', Z_dim(1) - 0.6, 'ind', 0)];
for l = 1:length(surf_p);
    [~, surf_p(l).ind] = min(abs(Z_dim - surf_p(l).z));
end;
lp = 0;
for l = 1:length(surf_p)
    for i = lb:xcs:rb
        lp = lp + 1;
        coordinates.surf(lp, 1) = XC(i);
        coordinates.surf(lp, 2) = YC(1);
        coordinates.surf(lp, 3) = Zr(surf_p(l).ind);
    end;
end;
Np = lp;
%* end

%** Abovebathymetry points calculating
%* begin
% There will be 4 levels of points as it is in the article.
bot_p = [struct('z', 0.1, 'ind', 0) ...
         struct('z', 0.2, 'ind', 0) ...
         struct('z', 0.4, 'ind', 0) ...
         struct('z', 0.6, 'ind', 0)];
for l = 1:length(bot_p);
    [~, bot_p(l).ind] = min(abs(Z_dim - bot_p(l).z));
end;
lp = 0;
for l = 1:length(bot_p)
    for i = lb:xcs:rb
        lp = lp + 1;
        coordinates.depth(lp, 1) = XC(i);
        coordinates.depth(lp, 2) = YC(1);
        coordinates.depth(lp, 3) = Zr(bot_p(l).ind);
    end;
end;
Np = Np + lp;
%* end

%** Onpycnocline points calculating
%* begin
lp = 0;
[~, indx] = min(abs(XC - h.Lg));
[~, indzh] = min(abs(-Zr + h.zpl));
[~, indzl] = min(abs(-Zr + h.zpr));
for i = lb:xcs:rb
    lp = lp + 1;
    coordinates.pyc(lp, 1) = XC(i);
    coordinates.pyc(lp, 2) = YC(1);
    if (i < indx)
        coordinates.pyc(lp, 3) = Zr(indzh);
    elseif (i >= indx)
        coordinates.pyc(lp, 3) = Zr(indzl);
    end
end;
for k=indzl+1:zcs:indzh
    lp = lp + 1;
    coordinates.pyc(lp, 1) = XC(indx);
    coordinates.pyc(lp, 2) = YC(1);
    coordinates.pyc(lp, 3) = Zr(k);
end;
Np = Np + lp;
%* end

%** Depict points
%* begin
if (h.deb ~= 0)
    figure(1);
    plot(coordinates.surf(:, 1), coordinates.surf(:, 3), 'r*');
    hold on;
    plot(coordinates.depth(:, 1), coordinates.depth(:, 3), 'r*');
    plot(coordinates.pyc(:, 1), coordinates.pyc(:, 3), 'r*');
    axis([0 h.Lx -h.H 0]);
    hold off;
end;
%* end

%** File structure:
%* See pkg/flt/README.flt for details
%* 1 - npart   A unique float identifier (1,2,3,...)
%* 2 - tstart  start date of integration of float (in s)
%*             Note: If tstart=-1 floats are integrated right from the
%*             beginning
%* 3 - xpart   x position of float (in units of XC)
%* 4 - ypart   y position of float (in units of YC)
%* 5 - kpart   actual vertical level of float (in units of vertical coordinates)
%* 6 - kfloat  target level of float (should be the same as kpart at
%*     the beginning)
%* 7 - iup     flag if the float
%*              - should profile   ( >  0 = return cycle (in s) to surface)
%*              - remain at depth  ( =  0 )
%*              - is a 3D float    ( = -1 ).
%*              - should be advected WITHOUT additional noise ( = -2 ).
%*                (This implies that the float is non-profiling)
%*              - is a mooring     ( = -3 ), i.e. the float is not advected
%* 8 - itop    time of float the surface (in s)
%* 9 - tend    end  date of integration of float (in s)
%*             Note: If tend=-1 floats are integrated till the end of
%*             the integration
%*
%* In addition the first line of the file contains a record with
%*      - the number of floats on that tile in the first record
%*      - the total number of floats in the sixth record
%* Note: At first-time initialization in a global file both fields
%* should be the same.

%** Trajectories parameters
%* begin
llength = 9;
%Np =;
tstart = -1.0;
%xpart =
%ypart =
%kpart = ;
%kfloat = kpart;
iup = -1.0;
itop = 0.0;
tend = -1.0;
%* end

%* First field !HEADER!
%*  1    2    3    4    5    6   7   8   9
%* [Np   0    0    0    0   Np   0   0   0]
%* Rest fields !FIELDS!
%*    1      2     3     4     5      6    7    8    9
%* [npart|tstart|xpart|ypart|kpart|kfloat|iup|itop|tend]

%** Output array making
%* begin
output = zeros(llength, Np + 1);
output(:, 1) = [Np 0 0 0 0 Np 0 0 0];
np = 1;
for i = 1:size(coordinates.surf, 1)
    output(:, np + 1) = [np  tstart coordinates.surf(i, 1) ...
                                    coordinates.surf(i, 2) ...
                                    coordinates.surf(i, 3) ...
                                    coordinates.surf(i, 3) ...
                         iup itop   tend];
    np = np + 1;
end;
for i = 1:size(coordinates.depth, 1)
    output(:, np + 1) = [np  tstart coordinates.depth(i, 1) ...
                                    coordinates.depth(i, 2) ...
                                    coordinates.depth(i, 3) ...
                                    coordinates.depth(i, 3) ...
                         iup itop   tend];
    np = np + 1;
end;
for i = 1:size(coordinates.pyc, 1)
    output(:, np + 1) = [np  tstart coordinates.pyc(i, 1) ...
                                    coordinates.pyc(i, 2) ...
                                    coordinates.pyc(i, 3) ...
                                    coordinates.pyc(i, 3) ...
                         iup itop   tend];
    np = np + 1;
end;
%* end

%** Save data
%* begin
fid=fopen(outfilename, 'w', ieee); fwrite(fid, output, prec); fclose(fid);
%* end

clear prec ieee outfilename depthfilename Lx Ly Lz Lg nx ny nz dx ...
    dy dy dz dz depth XC YC Np np Zr tstart kpart kfloat ...
    iup itop tend output fid i j k coordinates amplitude delta ind kpart_surf ...
    llength lp xcs zcs ans hpl hpr indx indzh indzl lb rb h CWD par_path ...
    Z_dim bot_p surf_p l;
