% Dimensionless horizontal velocity picture genaration.
clear CWD T X X_dim Z Z_dim dir_path filenames h nz outname par_path pic ...
    point rev t_b z_b  xTick yTick xTickMarks yTickMarks scale;

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
outname = 'reverse_jet.tiff';
% Area parameters reading
par_path = '../input/shares';
% The scale of the pictires
scale = 1.3;
CWD = pwd;
cd(par_path);
h = parameters;
cd(CWD);
%* end

%** Data reading and consolidating
%* begin
X = nc_varget(filenames{1}, 'X');
Z = nc_varget(filenames{1}, 'Z');
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
point = struct('x_p', 0.66, ...
               'x_i', 0,    ...
               'z_b', 0.5,  ...
               'z_i', 0);
[~, point.x_i] = min(abs(X_dim - point.x_p));
[~, point.z_i] = min(abs(Z_dim - point.z_b));
%* end

%** Dpeth of reversal jet finding
%* begin
rev = [];
for i = 1:length(filenames)
    u = nc_varget(filenames{i}, 'U')/h.c0;
    [~, tmp] = min(abs(u(:, point.z_i:1:nz, point.x_i)), [], 2);
    rev = [rev; tmp];
end;
rev = rev + point.z_i - 1;
rev = Z_dim(rev);
clear u i tmp;
%* end

%** Result drawing
%* begin
% 150% - [300,250,494,275]
% 200% - [300,250,674,360]
pic = figure('Position', [0,0,494*scale,275*scale]);
set(pic, 'PaperPositionMode', 'auto');
set(gca, 'FontSize', 16);
t_b = [85 115];
z_b = [0.0 0.45];
plot(T, rev, 'v');
grid on;
yTick = [0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45];
yTickMarks = {'0.0' '' '0.1' '' '0.2' '' '0.3' '' '0.4' ''};
xTick = [85 90 96 100 105 110 115];
xTickMarks = {'' '90' '' '100' '' '110' ''};
set(gca, 'ytick', yTick, 'YTickLabel', yTickMarks);
set(gca, 'xtick', xTick, 'XTickLabel', xTickMarks);
axis([t_b(1) t_b(2) z_b(1) z_b(2)]);
xlabel('tc_{0}/h_{2}','FontSize',18);
ylabel('z/h_{2}','FontSize',18);
print(pic, '-dtiff', '-r95', outname);
%* end

clear CWD T X X_dim Z Z_dim dir_path filenames h nz outname par_path pic ...
    point rev t_b z_b xTick yTick xTickMarks yTickMarks scale;
