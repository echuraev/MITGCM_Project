% Dimensionless interfacial displacement picture genaration.
clear dir_path filenames par_path CWD h pic X Z T Z_dim X_dim x_point ...
    z_point delta zindex xindex rho rho_disp eta_dim t_g t_b ...
    eta_b Tint eta_dim_int yTick yTickMarks i k t outname nz ...
    xTick xTickMarks scale;

%** Main variable definition
%* begin
dir_path = '../build/mnc_test_0001/';
filenames = {[dir_path 'state.0000000000.t001.nc']...
             [dir_path 'state.0000002824.t001.nc']};
outname = 'interfacial_dicplacement.tiff';
% Area parameters reading
par_path = './shares';
% The scale of the pictires
scale = 1.1;
CWD = pwd;
cd(par_path);
h = parameters;
cd(CWD);
%* end

%** Data reading and consolidation
%* begin
X = nc_varget(filenames{1}, 'X');
Z = nc_varget(filenames{1}, 'Z');
T = nc_varget(filenames{1}, 'T')*h.c0/h.h2r;
for i=2:length(filenames)
    T = cat(1, T, nc_varget(filenames{i}, 'T')*h.c0/h.h2r);
end;
%* end

%** Dimensionless variables
%* begin
X_dim = (X - h.Lg)/h.Lx;
x_point = 0.66;
Z_dim = (h.H + Z)/h.h2r;
%* end

%** Indexes of displacement point finding
%* begin
[~, xindex] = min(abs(X_dim - x_point));
%* end

%** Interfacial displacement finding
%* begin
eta_dim = zeros(length(T));
t_g = 1;
for i=1:length(filenames)
    rho = nc_varget(filenames{i}, 'Temp');
    for t=1:size(rho, 1)
        [~, zindex] = min(abs(rho(t, :, xindex)));
        % zindex has to be a bit nearer to the real one.
        eta_dim(t_g) = Z_dim(zindex + 2);
        t_g = t_g + 1;
    end;
    clear rho;
end;
%* end

%** Result drawing
%* begin
% 150% - [300,250,515,335]
% 200% - [300,250,645,440]
pic = figure('Position', [0,0,692*scale,447*scale]);
set(pic, 'PaperPositionMode', 'auto');
set(gca, 'FontSize', 16);
t_b = [50 100];
eta_b = [2.4 4.4];
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
axis([t_b(1) t_b(2) eta_b(1) eta_b(2)]);
xlabel('tc_{0}/h_{2}','FontSize',18);
ylabel('\eta/h_{2}','FontSize',18);
print(pic, '-dtiff', '-r95', outname);
%* end

clear dir_path filenames par_path CWD h pic X Z T Z_dim X_dim x_point ...
    z_point delta zindex xindex rho rho_disp eta_dim t_g t_b ...
    eta_b Tint eta_dim_int yTick yTickMarks i k t outname nz ...
    xTick xTickMarks scale;
