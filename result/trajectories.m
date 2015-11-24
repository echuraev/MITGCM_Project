% Depict the groups of particles trajectories.
clear CWD Np_bot Np_pyc Np_surf colors data dir_path filenames flt grid_x ...
    grid_z h header i np outname par_path points traj_dir x_b z_b ...
    bot_points pyc_points surf_points cl1 cl2 l leg pic pl st en scale ind;

%** Input parameters
%* begin
outname = 'trajectories.tiff';
% Area parameters reading
par_path = '../input/shares';
% The scale of the pictires
scale = 1.1;
ttime = 57; % s
CWD = pwd;
cd(par_path);
h = parameters;
cd(CWD);
%------------------------
traj_dir = '../build/';
%* end

%** Data reading
%* begin
[flt,~,~] = read_flt_traj(traj_dir);
%* end

%** Depict trajectories only up to 57 s. moments. Cut uninteresting points
%** for now.
%* begin
[~, ind] = min(abs(flt(1,1).time - ttime));
for l=1:length(flt)
    flt(1, l).time(ind:end) = [];
    flt(1, l).x(ind:end) = [];
    flt(1, l).y(ind:end) = [];
    flt(1, l).z(ind:end) = [];
    flt(1, l).i(ind:end) = [];
    flt(1, l).j(ind:end) = [];
    flt(1, l).k(ind:end) = [];
    flt(1, l).p(ind:end) = [];
    flt(1, l).u(ind:end) = [];
    flt(1, l).v(ind:end) = [];
    flt(1, l).t(ind:end) = [];
    flt(1, l).s(ind:end) = [];
end;
%* end

%** Depicting points initialization
%* begin
surf_points = [
    struct('name', '5.2', 'arr', [], 'st', 3, 'gap', 20, 'Np', 100, 'color', '-b'); ... % surface
    struct('name', '5.1', 'arr', [], 'st', 5, 'gap', 20, 'Np', 100, 'color', '-m'); ... % surface
    struct('name', '4.9', 'arr', [], 'st', 7, 'gap', 20, 'Np', 100, 'color', '-y'); ... % surface
    struct('name', '4.7', 'arr', [], 'st', 9, 'gap', 20, 'Np', 100, 'color', '-k')
];    % surface
               
bot_points = [
    struct('name', '0.1', 'arr', [], 'st', 5, 'gap', 10, 'Np', 100, 'color', '-b'); ... % bottom
    struct('name', '0.2', 'arr', [], 'st', 5, 'gap', 10, 'Np', 100, 'color', '-m'); ... % bottom
    struct('name', '0.4', 'arr', [], 'st', 5, 'gap', 10, 'Np', 100, 'color', '-y'); ... % bottom
    struct('name', '0.6', 'arr', [], 'st', 5, 'gap', 10, 'Np', 100, 'color', '-k')
]; % bottom
pyc_points = [
    struct('name', 'on-pycnocline', 'arr', [], 'st', 3, 'gap', 7, 'Np', 142, 'color', '-b')
];    % pycline
%* end

%** Start points arrays building
%* begin
st = 1;
en = 0;
for l=1:length(surf_points)
    en = en + surf_points(l).Np;
    surf_points(l).arr = flt(st:en);
    st = st + surf_points(l).Np;
end;
for l=1:length(bot_points)
    en = en + bot_points(l).Np;
    bot_points(l).arr = flt(st:en);
    st = st + bot_points(l).Np;
end;
for l=1:length(pyc_points)
    en = en + pyc_points(l).Np;
    pyc_points(l).arr = flt(st:en);
    st = st + pyc_points(l).Np;
end;
%* end

%** Depict resulting arrayes
%* begin
% Undersurface points depicting
pic = figure('Position', [0,0,1280*scale,250*scale], 'PaperPositionMode','auto');
set(gca, 'FontSize', 16);
z_b = [2.0 5.4];
x_b = [-0.1 0.9];
leg.lhs = zeros(1, length(surf_points) + 2);
leg.names = cell(1, length(surf_points) + 2);
hold on;
for np = 1:length(surf_points)
    for i = surf_points(np).st:surf_points(np).gap:length(surf_points(np).arr)
        % plot curves
        pl = plot((surf_points(np).arr(i).x - h.Lg)/h.Lx, ...
                  (surf_points(np).arr(i).z + h.H)/h.h2r, ...
                  surf_points(np).color);
        % Mark bagining and ending of cureves
        cl1 = plot((surf_points(np).arr(i).x(1) - h.Lg)/h.Lx, (surf_points(np).arr(i).z(1) + h.H)/h.h2r, 'g*');
        cl2 = plot((surf_points(np).arr(i).x(end) - h.Lg)/h.Lx, (surf_points(np).arr(i).z(end) + h.H)/h.h2r, 'r*');
    end;
    % Add legend for the color;
    leg.lhs(np) = pl;
    leg.names{np} = surf_points(np).name;
end;
leg.lhs(end-1) = cl1;
leg.names{end-1} = 'initial position';
leg.lhs(end) = cl2;
leg.names{end} = 'final position';
legend(leg.lhs(:), leg.names{:}, 'Orientation', 'horizontal', ...
       'Location', 'SouthWest');
hold off;
grid on;
axis([x_b(1) x_b(2) z_b(1) z_b(2)]);
xlabel('x/L_{x}', 'FontSize', 18);
ylabel('z/h_{2}', 'FontSize', 18);
print(pic, '-dtiff', '-r95', ['surf_' outname]);

% Abovebottom points depicting
pic = figure('Position', [0,0,1280*scale,250*scale], 'PaperPositionMode','auto');
set(gca, 'FontSize', 16);
z_b = [0.0 1.1];
x_b = [-0.1 0.9];
leg.lhs = zeros(1, length(bot_points) + 2);
leg.names = cell(1, length(bot_points) + 2);
hold on;
for np = 1:length(bot_points)
    for i = bot_points(np).st:bot_points(np).gap:length(bot_points(np).arr)
        % plot curves
        pl = plot((bot_points(np).arr(i).x - h.Lg)/h.Lx, ...
                  (bot_points(np).arr(i).z + h.H)/h.h2r, ...
                  bot_points(np).color);
        % Mark bagining and ending of cureves
        cl1 = plot((bot_points(np).arr(i).x(1) - h.Lg)/h.Lx, (bot_points(np).arr(i).z(1) + h.H)/h.h2r, 'g*');
        cl2 = plot((bot_points(np).arr(i).x(end) - h.Lg)/h.Lx, (bot_points(np).arr(i).z(end) + h.H)/h.h2r, 'r*');
    end;
    % Add legend for the color;
    leg.lhs(np) = pl;
    leg.names{np} = bot_points(np).name;
end;
leg.lhs(end-1) = cl1;
leg.names{end-1} = 'initial position';
leg.lhs(end) = cl2;
leg.names{end} = 'final position';
legend(leg.lhs(:), leg.names{:}, 'Orientation', 'horizontal', ...
       'Location', 'North');
hold off;
grid on;
axis([x_b(1) x_b(2) z_b(1) z_b(2)]);
xlabel('x/L_{x}', 'FontSize', 18);
ylabel('z/h_{2}', 'FontSize', 18);
print(pic, '-dtiff', '-r95', ['bot_' outname]);

% Onpycnocline points depicting
pic = figure('Position', [0,0,1280*scale,250*scale], 'PaperPositionMode','auto');
set(gca, 'FontSize', 16);
z_b = [1.0 5.0];
x_b = [-0.1 0.9];
leg.lhs = zeros(1, length(pyc_points) + 2);
leg.names = cell(1, length(pyc_points) + 2);
hold on;
for np = 1:length(pyc_points)
    for i = pyc_points(np).st:pyc_points(np).gap:length(pyc_points(np).arr)
        % plot curves
        pl = plot((pyc_points(np).arr(i).x - h.Lg)/h.Lx, ...
                  (pyc_points(np).arr(i).z + h.H)/h.h2r, ...
                  pyc_points(np).color);
        % Mark bagining and ending of cureves
        cl1 = plot((pyc_points(np).arr(i).x(1) - h.Lg)/h.Lx, (pyc_points(np).arr(i).z(1) + h.H)/h.h2r, 'g*');
        cl2 = plot((pyc_points(np).arr(i).x(end) - h.Lg)/h.Lx, (pyc_points(np).arr(i).z(end) + h.H)/h.h2r, 'r*');
    end;
    % Add legend for the color;
    leg.lhs(np) = pl;
    leg.names{np} = pyc_points(np).name;
end;
leg.lhs(end-1) = cl1;
leg.names{end-1} = 'initial position';
leg.lhs(end) = cl2;
leg.names{end} = 'final position';
legend(leg.lhs(:), leg.names{:}, 'Orientation', 'horizontal', ...
       'Location', 'South');
hold off;
grid on;
axis([x_b(1) x_b(2) z_b(1) z_b(2)]);
xlabel('x/L_{x}', 'FontSize', 18);
ylabel('z/h_{2}', 'FontSize', 18);
print(pic, '-dtiff', '-r95', ['pyc_' outname]);
%* end

clear CWD Np_bot Np_pyc Np_surf colors data dir_path filenames flt grid_x ...
    grid_z h header i np outname par_path points traj_dir x_b z_b ...
    bot_points pyc_points surf_points cl1 cl2 l leg pic pl st en scale ind;
