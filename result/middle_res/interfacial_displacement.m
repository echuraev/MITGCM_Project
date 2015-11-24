% Dimensionless interfacial displacement picture genaration.
clear filenames Hr H Lx Ly Lg Lw r2 r1 h2r h1r V a c_max dh g nu i ReW T X ...
    X_dim Z c0 eta_dim eta_diml eta_dimh xindex k nx nz r_m t t_diml rho ...
    rho_disp t_dimr t_g x_point h Z_dim z_point zindex dir_name dir_prefix ...
    exp_tar Tint eta_dim_int xTick xTickMarks yTick yTickMarks Hl dH dSurf ...
    delta dg dr h1 h1l h2 h2add h2l zpl zpr;

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
%c0 = 0.112;
c_max = sqrt(g*(h1r + h2r)*(r1 - r2)/((sqrt(r1) + sqrt(r2))^2));
% Internal soliton amplitude
a = (h1r*sqrt(r2) - h2r*sqrt(r1))/(sqrt(r1) + sqrt(r2));
% Wave Reynolds number
ReW = c0*H/nu; % ~2.6*10e4
%* end

%** Data reading and consolidation
%* begin
X = nc_varget(filenames{1}, 'X');
Z = nc_varget(filenames{1}, 'Z');
nx = length(X);
nz = length(Z);
%rho = r_m - nc_varget(filenames{1}, 'Temp');
T = nc_varget(filenames{1}, 'T')*c0/h2r;
for i=2:length(filenames)
%    rho = cat(1, rho, (r_m - nc_varget(filenames{i}, 'Temp')));
    T = cat(1, T, nc_varget(filenames{i}, 'T')*c0/h2r);
end;
%* end

%** Dimensionless variables
%* begin
X_dim = (X - Lg)/Lx;
x_point = 0.66;
Z_dim = (H + Z)/h2r;
% zpr - depth of interface.
z_point = (H + zpr)/h2r; % 4.2
%* end

%** Indexes of displacement point finding
%* begin
[delta xindex] = min(abs(X_dim - x_point));
[delta zindex] = min(abs(Z_dim - z_point));
zindex = zindex + 1;
%* end

%** Interfacial displacement finding
%* begin
rho = r_m - nc_varget(filenames{1}, 'Temp');
rho_disp = rho(1, zindex, xindex);
eta_dim(length(T)) = 0.0;
t_g = 1;
for i=2:length(filenames)
    for t=1:size(rho, 1)
        for k=1:nz
            if (rho(t, k, xindex) >= rho_disp)
                eta_dim(t_g) = Z_dim(k);
                break;
            end;
        end;
        t_g = t_g + 1;
    end;
    clear rho;
    rho = r_m - nc_varget(filenames{i}, 'Temp');
end;
%* end

%** Result drawing
%* begin
h = figure('Position', [300,250,515,335]);
set(h,'PaperPositionMode','auto');
t_diml = 50;
t_dimr = 100;
eta_dimh = 4.4;
eta_diml = 2.4;
% May be interpolaton necessary.
Tint = T(1):1.0:T(size(T), 1);
eta_dim_int = interp1(T, eta_dim, Tint);
plot(Tint, eta_dim_int, 'LineWidth', 2);
grid on;
yTick = 2.4:0.2:4.4;
xTick = 50:5:100;
yTickMarks = {};
xTickMarks = {};
for i=1:size(yTick, 2)
    yTickMarks(i) = {num2str(yTick(i))};
end;
for i=1:size(xTick, 2)
    xTickMarks(i) = {num2str(xTick(i))};
end;
set(gca, 'ytick', yTick, 'YTickLabel', yTickMarks);
set(gca, 'xtick', xTick, 'XTickLabel', xTickMarks);
axis([t_diml t_dimr eta_diml eta_dimh]);
xlabel('tc_{0}/h_{2}','FontSize',18);
ylabel('\eta/h_{2}','FontSize',18);
print(h, '-dtiff', '-r95', '02_interfacial_displacement.tiff');
%* end

clear filenames Hr H Lx Ly Lg Lw r2 r1 h2r h1r V a c_max dh g nu i ReW T X ...
    X_dim Z c0 eta_dim eta_diml eta_dimh xindex k nx nz r_m t t_diml rho ...
    rho_disp t_dimr t_g x_point h Z_dim z_point zindex dir_name dir_prefix ...
    exp_tar Tint eta_dim_int xTick xTickMarks yTick yTickMarks Hl dH dSurf ...
    delta dg dr h1 h1l h2 h2add h2l zpl zpr;
