% Dimensionless horizontal velocity picture genaration.
clear dir_path filenames par_path CWD h pic X Z T Z_dim X_dim x_point ...
    nz_points z_points delta u_dim t_g t_b u_b ...
    Tint u_dim_int yTick yTickMarks i p t outname scale;

%** Main variable definition
%* begin
dir_path = '../build/mnc_test_0001/';
filenames = {[dir_path 'state.0000000000.t001.nc']...
             [dir_path 'state.0000002824.t001.nc']};
outname = 'horizontal_velocity.tiff';
% Area parameters reading
par_path = './shares';
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
T = nc_varget(filenames{1}, 'T')*h.c0/h.h2r;
for i=2:length(filenames)
    T = cat(1, T, nc_varget(filenames{i}, 'T')*h.c0/h.h2r);
end;
%* end

%** Dimensionless variables
%* begin
Z_dim = (h.H + Z)/h.h2r;
X_dim = (X - h.Lg)/h.Lx;
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
u_dim(nz_points, length(T)) = 0.0;
t_g = 1;
for i=1:length(filenames)
    u = nc_varget(filenames{i}, 'U');
    for p=1:nz_points
        for t=1:size(u, 1)
            u_dim(p, t_g + t - 1) = u(t, z_points(p, 2), x_point(2))/h.c0;
        end;
    end;
    t_g = t_g + size(u, 1);
    clear u
end;
%* end

%** Result drawing
%* begin
% 150% - [300,250,515,335]
% 200% - [300,250,684,439]
pic = figure('Position', [0,0,684*scale,439*scale]);
set(pic, 'PaperPositionMode', 'auto');
set(gca, 'FontSize', 16);
t_b = [50 110];
u_b = [-0.8 0.3];
% May be interpolaton necessary.
Tint = T(1):0.7:T(size(T,1));
for i=1:size(u_dim, 1)
    u_dim_int(i, :) = interp1(T, u_dim(i, :), Tint);
end;
plot(Tint, u_dim_int(1, :), '+');
grid on;
yTick = [-0.7 -0.5 -0.3 -0.1 0.1 0.3];
yTickMarks = {'-0.7' '-0.5' '-0.3' '-0.1' '0.1' '0.3'}; 
hold on;
plot(Tint, u_dim_int(2, :), 'x');
plot(Tint, u_dim_int(3, :), '.');
plot(Tint, u_dim_int(4, :), 'o');
set(gca, 'ytick', yTick, 'YTickLabel', yTickMarks);
legend('0.1', '0.2', '0.4', '0.6', 'Location', 'SouthEast');
hold off;
axis([t_b(1) t_b(2) u_b(1) u_b(2)]);
xlabel('tc_{0}/h_{2}','FontSize',18);
ylabel('u/c_{0}','FontSize',18);
print(pic, '-dtiff', '-r95', outname);
%* end

clear dir_path filenames par_path CWD h pic X Z T Z_dim X_dim x_point ...
    nz_points z_points delta u_dim t_g t_b u_b ...
    Tint u_dim_int yTick yTickMarks i p t outname scale;
