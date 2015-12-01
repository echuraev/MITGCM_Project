clear p;

p.build_path = '../build/';
p.outname = 'U.gif';

p.U = squeeze(rdmds ([p.build_path 'U'], NaN));
p.W = squeeze(rdmds ([p.build_path 'W'], NaN));
p.X = squeeze(rdmds ([p.build_path 'XG']));
p.Y = squeeze(rdmds ([p.build_path 'RC']));
%p.Depth = squeeze(rdmds ([p.build_path 'Depth']));
[p.X1, p.Y1] = meshgrid (p.X, p.Y);
p.U(p.U == 0) = NaN;
p.W(p.W == 0) = NaN;
p.Xvector = p.X(1:50:size(p.X,1),1);
p.Yvector = p.Y(1:20:size(p.Y,1),1);
[p.Xvector, p.Yvector] = meshgrid (p.Xvector, p.Yvector);
p.Uvector = p.U(1:50:size(p.U,1),1:20:size(p.U,2),:);
p.Wvector = p.W(1:50:size(p.W,1),1:20:size(p.W,2),:);
%p.h = figure(1);
p.h = figure('Position', [100 100 1200 300]);
p.minU = min(min(min(sqrt(p.U.^2 + p.W.^2))));
p.maxU = max(max(max(sqrt(p.U.^2 + p.W.^2))));

% уменьшить длинну вектора, пропорциональные размеры, стрелочки
%set
for i = 1:size(p.U, 3)    
    pcolor (p.X1', p.Y1', sqrt(p.U(:,:,i).^2 + p.W(:,:,i).^2));
    caxis ([p.minU p.maxU]);
    hold on;
    quiver(p.Xvector',p.Yvector',p.Uvector(:,:,i),p.Wvector(:,:,i), 0, 'color',[0 0 0]);
    %plot(p.X, -p.Depth);
    shading interp;
    p.cb = colorbar;
    xlabel ('x, m'); %help
    ylabel ('y, m'); %help
    title ('U');
    title (p.cb, 'm/s');
    hold off;
    drawnow;
    p.f = getframe(p.h);
    p.im = frame2im(p.f);
    [p.gif, p.map] = rgb2ind(p.im, 256);
    if (i == 1)
        imwrite(p.gif, p.map, p.outname, 'DelayTime', 0.3, 'LoopCount', inf);
    else
        imwrite(p.gif, p.map, p.outname, 'WriteMode', 'append');
    end;
    
    
end

clear p;
