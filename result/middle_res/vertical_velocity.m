% Dimensionless vertical velocity picture genaration.
clear filenames Hr H Lx Ly Lg Lw r2 r1 h2r h1r V h1 h2 g nu i ReW T X X_dim ...
    Z c0 w_dim w_diml w_dimh index k nx nz r_m t t_diml t_dimr t_g x_point ...
    Z_dim nz_points z_points p dir_name dir_prefix exp_tar yTickMarks yTick ...
    w_dim_int Tint Hl dH dSurf delta dg dh dr h h1l h2add h2l zpl zpr;

%** Main variable definition
%* begin
dir_prefix='../build/';
dir_name='mnc_test_0001';
exp_tar='';
filenames={[dir_prefix dir_name exp_tar '/state.0000000000.t001.nc']...
           [dir_prefix dir_name exp_tar '/state.0000004236.t001.nc']};

%** Area parameters calculating using algorithm in section above.
%* begin
% Size of X domain
Lx=6.4;
% Size of Y domain
Ly=0.4;
%* Length of the domain to the left of the gate
Lg = 0.6; % meters
%* Width of the domain
Lw = 0.4; % meters
% Density of the layers
r1 = 1047.0;
r2 = 1022.0;
% Volume of added water
V=38;
% Initial water levels
h1 = 0.2;
h2 = 0.05;
% Depth of the added liquid.
h2add = (V/r2)/(Lg*Lw);
% h2add  = 0.158; - In article
% h2add  = 0.1549; - We have
% Target h2 levels calculating.
h2r = h2;
h2l = h2 + h2add;
%h2r = 0.05; % meters - In article
%h2r = 0.05; % meters - We have
%h2l = 0.21; % meters - In article
%h2l = 0.2049; % meters - We have (It gives less amplitude of a soliton)
% Delta of r1 levels calculating.
dH = r2*(h2l - h2)/r1;
% r1 level to the left and to the right of gate calculating.
h1l = h1 - dH*(1 - Lg/Lx);
h1r = h1l + dH;
% The target depth of water to the laft and to the rigtht of gate
% calculating.
Hl = h2l + h1l;
Hr = h2r + h1r; % or Hr = h2 + h1r;
H = Hl;
% Common level is H, but surface to the left of the gate is 
% increased on (Hl - Hr). And depth of the left pycnocline depth will be
% decreased.
dSurf = Hl - Hr;
%h2l = h2l - dSurf;
% Rest of parameters setting
% Pycnocline thickness
dh=0.015; % meters - In articlce
% Gate thickness
% But the real thickness is 0.005.
dg=0.01; % meters % may be 0.01 meters because mixing
% Pycnoclines depth
zpl=-h2l;
zpr=-h2r;
% Middle and delta of density
r_m=(r1 + r2)/2;
dr=(r1 - r2);
% Total acceleration and kinematic viscosity
g = 9.81;
nu = 1.e-6;
%* end

%** Additional constants calculations.
%* begin
% internal long wave velocity
c0 = sqrt(g*((r1 - r2)/(r1))*((h1r*h2r)/(h1r + h2r)));
% Wave Reynolds number
ReW = c0*H/nu; % ~2.6*10e4
%* end

%** Data reading and consolidation
%* begin
X = nc_varget(filenames{1}, 'X');
Z = nc_varget(filenames{1}, 'Z');
nx = length(X);
nz = length(Z);
T = nc_varget(filenames{1}, 'T')*c0/h2r;
for i=2:length(filenames)
    T = cat(1, T, nc_varget(filenames{i}, 'T')*c0/h2r);
end;
%* end

%** Dimensionless variables
%* begin
Z_dim = (H + Z)/h2r;
X_dim = (X - Lg)/Lx;
x_point = [0.66 0];
nz_points = 4;
z_points = [0.1 0; 0.2 0; 0.4 0; 0.6 0];
%* end

%** Indexes of veloctity points finding.
%* begin
[delta x_point(2)] = min(abs(X_dim - x_point(1)));
for p=1:nz_points
    [delta z_points(p, 2)] = min(abs(Z_dim - z_points(p, 1)));
end;
%* end

%** Interfacial displacement finding
%* begin
w_dim(nz_points, length(T)) = 0.0;
t_g = 1;
for i=1:length(filenames)
    w = nc_varget(filenames{i}, 'W');
    for p=1:nz_points
        for t=1:size(w, 1)
            w_dim(p, t_g + t - 1) = w(t, z_points(p, 2), x_point(2))/c0;
        end;
    end;
    t_g = t_g + size(w, 1);
    clear w
end;
%* end

%** Result drawing
%* begin
h = figure('Position', [300,250,515,335]);
set(h,'PaperPositionMode','auto');
t_diml = 50;
t_dimr = 110;
w_dimh = 0.05;
w_diml = -0.05;
% May be interpolaton necessary.
Tint = T(1):0.7:T(size(T,1));
for i=1:size(w_dim, 1)
    w_dim_int(i, :) = interp1(T, w_dim(i, :), Tint);
end;
plot(Tint, w_dim_int(1, :), '+');
grid on;
yTick = [-0.04 -0.02 0.0 0.02 0.04];
yTickMarks = {'-0.04' '-0.02' '0.0' '0.02' '0.04'}; 
hold on;
plot(Tint, w_dim_int(2, :), 'x');
plot(Tint, w_dim_int(3, :), '.');
plot(Tint, w_dim_int(4, :), 'o');
set(gca, 'ytick', yTick, 'YTickLabel', yTickMarks);
legend('0.1', '0.2', '0.4', '0.6', 'Location', 'SouthEast');
hold off;
axis([t_diml t_dimr w_diml w_dimh]);
xlabel('tc_{0}/h_{2}','FontSize',18);
ylabel('w/c_{0}','FontSize',18);
print(h, '-dtiff', '-r95', '03_vertical_velocity.tiff');
%* end

clear filenames Hr H Lx Ly Lg Lw r2 r1 h2r h1r V h1 h2 g nu i ReW T X X_dim ...
    Z c0 w_dim w_diml w_dimh index k nx nz r_m t t_diml t_dimr t_g x_point ...
    Z_dim nz_points z_points p dir_name dir_prefix exp_tar yTickMarks yTick ...
    w_dim_int Tint Hl dH dSurf delta dg dh dr h h1l h2add h2l zpl zpr;
