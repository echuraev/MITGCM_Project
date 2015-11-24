% Dimensionless density depicting.
clear CWD T X X_dim Z Z_dim dir_path f_index filenames grid_x grid_z h i ...
    nx nz outname par_path pic rho t_index t_point x_b x_point z_b ...
    d_c map max_c min_c cb scale;

%** Main variable definition
%* begin
dir_path = '../build/mnc_test_0001/';
filenames = {[dir_path 'state.0000000000.t001.nc']...
             [dir_path 'state.0000001408.t001.nc']...
             [dir_path 'state.0000002816.t001.nc']...
             [dir_path 'state.0000004224.t001.nc']};%...
%             [dir_path 'state.0000005632.t001.nc']...
%             [dir_path 'state.0000007040.t001.nc']...
%             [dir_path 'state.0000008448.t001.nc']...
%             [dir_path 'state.0000009856.t001.nc']...
%             [dir_path 'state.0000011265.t001.nc']...
%             [dir_path 'state.0000012672.t001.nc']};
outname = 'density.tiff';
% Area parameters reading
par_path = '../input/shares';
% The scale of the pictires
scale = 1.1;
CWD = pwd;
cd(par_path);
h = parameters;
cd(CWD);
%* end

%** Data reading and consolidating
%* begin
X = nc_varget(filenames{1}, 'X');
Z = nc_varget(filenames{1}, 'Z');
nx = length(X);
nz = length(Z);
T = nc_varget(filenames{1}, 'T')*h.c0/h.h2r;
for i=2:length(filenames)
    T = cat(1, T, nc_varget(filenames{i}, 'T')*h.c0/h.h2r);
end;
%* end

%** Dimensionless variables
%* begin
Z_dim = (h.H + Z)/h.h2r;
X_dim = (X - h.Lg)/h.Lx;
x_point = 0.66;
t_point = 76.6;
%* end

%** Density at the interesting moment finding.
%* begin
% Find file with the interesting time point.
for f = 1:length(filenames)
    tmp = nc_varget(filenames{f}, 'T')*h.c0/h.h2r;
    if (t_point >= tmp(1) && t_point <= tmp(size(tmp, 1)))
        [delta t_index] = min(abs(tmp - t_point));
        f_index = f;
        break;
    end;
end;
clear tmp f delta;
% Read the interesting density.
rho = nc_varget(filenames{f_index}, 'Temp');
rho = squeeze((h.r_m - rho(t_index, :, :))/h.r0);
%* end

%** Results drawing
%* begin
z_b = [0.0 5.1];
x_b = [0.45 0.85];%[0.0 0.9];
[grid_x grid_z] = meshgrid(X_dim, Z_dim);
% 200% - [300, 250, 710, 435]
% 150% - [300, 250, 550, 325]
pic = figure('Position', [0,0,710*scale,435*scale]);%[156,311,1081,297]);
set(pic, 'PaperPositionMode', 'auto');
set(gca, 'FontSize', 16);
pcolor(grid_x, grid_z, rho);
shading interp;
% Custom colormap creating;
min_c = 0.09; max_c = 0.45; d_c = (max_c - min_c)/63;
map = zeros(64, 3);
map(:, 1) = min_c:d_c:max_c;
map(:, 2) = min_c:d_c:max_c;
map(:, 3) = min_c:d_c:max_c;
colormap(map);
cb = colorbar;
set(cb, 'FontSize', 16);
hold on;
contour(grid_x, grid_z, rho, 2, '-w');
hold off;
xlabel('x/L_{x}', 'FontSize', 18);
ylabel('z/h_{2}', 'FontSize', 18);
axis([x_b(1) x_b(2) z_b(1) z_b(2)]);
print(pic, '-dtiff', '-r95', outname);
%* end

clear CWD T X X_dim Z Z_dim dir_path f_index filenames grid_x grid_z h i ...
    nx nz outname par_path pic rho t_index t_point x_b x_point z_b ...
    d_c map max_c min_c cb scale;
