clear p;

p.build_path = '../build/';
p.outname = 'T.gif';

p.T = squeeze(rdmds ([p.build_path 'T'], NaN));
p.X = squeeze(rdmds ([p.build_path 'XG']));
p.Y = squeeze(rdmds ([p.build_path 'RC']));
p.Depth = squeeze(rdmds ([p.build_path 'Depth']));
[p.X1, p.Y1] = meshgrid (p.X, p.Y);
p.T(p.T == 0) = NaN;

p.h = figure(1);
p.minT = min(min(min(p.T)));
p.maxT = max(max(max(p.T)));
%set
for i = 1:size(p.T, 3)    
    pcolor (p.X1', p.Y1', p.T(:,:,i));
    caxis ([p.minT p.maxT]);
    hold on;
    plot(p.X, -p.Depth);
    shading interp;
    p.cb = colorbar;
    xlabel ('x, m'); %help
    ylabel ('y, m'); %help
    title ('T');
    title (p.cb, 'C');
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
