% Dimensionless horizontal velocity picture genaration.
clear CWD X X_dim Z Z_dim dir_path filenames h i nz outname par_path pic ...
    points u_b xTick xTickMarks yTick yTickMarks z_b scale;

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
outname = 'bot_velocity.tiff';
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
nz = length(Z);
%* end

%** Dimensionless variables
%* begin
Z_dim = (h.H + Z)/h.h2r;
X_dim = (X - h.Lg)/h.Lx;
% Array of structures which store info about particular time-point.
points = [struct('t_p',  88.28, ... % 88.28 - In article
                 't_i',  0, ...
                 'x_p',  0.66, ...
                 'x_i',  0, ...
                 'file', 0, ...
                 'u',    zeros(nz, 1)) ...
          struct('t_p',  104.8, ... % 104.8 - In article
                 't_i',  0, ...
                 'x_p',  0.66, ...
                 'x_i',  0, ...
                 'file', 0, ...
                 'u',    zeros(nz, 1))];
%* end

%** Indexes of the interesting points finding.
%* begin
% Find files with the interesting time points.
for f = 1:length(filenames)
    % T is not repeated in every file.
    tmp = nc_varget(filenames{f}, 'T')*h.c0/h.h2r;
    for i = 1:size(points, 2)
        if (points(i).t_p >= tmp(1) && points(i).t_p <= tmp(size(tmp, 1)))
            [~, points(i).t_i] = min(abs(tmp - points(i).t_p));
            points(i).file = f;
        end;
    end;
end;
for i = 1:size(points, 2)
    [~, points(i).x_i] = min(abs(X_dim - points(i).x_p));
end;
clear tmp f i;
%* end

%** Velocity points finding
%* begin
for i = 1:size(points, 2)
    u = nc_varget(filenames{points(i).file}, 'U');
    points(i).u = squeeze(u(points(i).t_i, :, points(i).x_i)/h.c0);
end;
clear f u;
%* end

%** Result drawing
%* begin
% 150% - [300,250,494,268]
% 200% - [300,250,684,439]
pic = figure('Position', [0,0,661*scale,360*scale]);
set(pic, 'PaperPositionMode', 'auto');
set(gca, 'FontSize', 16);
z_b = [0 1];
u_b = [-0.20 0.15];
plot(points(1).u, Z_dim, '+');
hold on;
plot(points(2).u, Z_dim, 'o');
hold off;
grid on;
yTick = [0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
yTickMarks = {'0.0' '' '0.2' '' '0.4' '' '0.6' '' '0.8' '' '1.0'};
xTick = [-0.20 -0.15 -0.10 -0.05 0.0 0.05 0.10 0.15];
xTickMarks = {'' '-0.15' '' '-0.05' '' '0.05' '' '0.15'}; 
set(gca, 'ytick', yTick, 'YTickLabel', yTickMarks);
set(gca, 'xtick', xTick, 'XTickLabel', xTickMarks);
legend(['t*=' num2str(points(1).t_p)], ...
       ['t=' num2str(points(2).t_p)], 'Location', 'NorthEast');
axis([u_b(1) u_b(2) z_b(1) z_b(2)]);
xlabel('u/c_{0}','FontSize',18);
ylabel('z/h_{2}','FontSize',18);
print(pic, '-dtiff', '-r95', outname);
%* end

clear CWD X X_dim Z Z_dim dir_path filenames h i nz outname par_path pic ...
    points u_b xTick xTickMarks yTick yTickMarks z_b scale;
