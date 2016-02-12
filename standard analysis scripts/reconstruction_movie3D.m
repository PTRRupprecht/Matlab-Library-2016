
timetracesX_X_raw{pp} ROI_mapXX{pp}






%% plotting the yz-scanning trajectory

tz = linspace(98.5+15,-150-15,2312);
% tz = tz(end:-1:1);

ty = (-128:135)*0.7258;
ty2 = (-128:119)*0.7258;
ty = [ty,ty2,ty,ty2,ty,ty2,ty,ty2,ty].*(1-0.1./100.*tz);


figure, hold on;
for k = 1:9
    plot(ty(((k-1)*256+1):k*256),tz(((k-1)*256+1):k*256));
end


counter = 1;
for pp = 1:5
    for k = 1:max(ROI_mapXX{pp}(:))
        [x,y] = find(ROI_mapXX{pp} == k);
        if ~isempty(x)
            yy = mean(x); xx = mean(y);
            z_index = round(yy+512*(pp-1));
            coo_x(counter) = (xx-256)*0.7258*(1-0.1./100.*tz(z_index));
            if yy > 257
                yy = yy - 257;
            end
            coo_y(counter) = (yy - 128)*0.7258*(1-0.3./100.*tz(z_index));
            coo_z(counter) = -tz(z_index);
%             coordinates{counter} = [pos_x,pos_y,pos_z];
            counter = counter + 1;
            
            tt = timetracesX_X_raw{pp}(:,k);
%             tt = tt + max(min(tt),0) + 5;
            F0 = max(0,min(smooth(tt,500)));
            tt = (tt - F0)/(F0+30);
            timetrace(counter,:) = tt;
        end
    end
end
figure(2), imagesc(smoothed_timetracesX_cooked,[-1 20])
smoothed_timetracesX_cooked = conv2(timetrace,fspecial('gaussian',[1,9],4.5),'same');



%% pre 1
cmap = [linspace(0.66,1,800);linspace(0.66,0,800);linspace(0.66,0,800)]';
% smoothed_timetracesX_cooked = conv2(timetracesX_cooked,fspecial('gaussian',[17,1],7),'same');
[x,y,z] = sphere(10);
for timetT = 1:100
    G = figure('Position',[ 110 110 860 660],'Color',[1 1 1]);
    surf(x*0.1+400, y*0.1+400, z*0.1+400,ones(size(z))*600); hold on;
    surf(x*0.1+400, y*0.1+400, z*0.1+400,ones(size(z))*-100); hold on;
    for i = 1:numel(coo_x)
        if ~isnan(coo_x(i))
            r = 2;
            surf(x*r+coo_x(i), y*r+coo_y(i), z*r+coo_z(i),ones(size(z))*600); hold on;
        end
    end
    colormap(cmap);
    caxis([-100 600])
    colorbar;
    shading(gca,'interp')
    hold off;
    axis([min(coo_x) max(coo_x) min(coo_y) max(coo_y) min(coo_z) max(coo_z)])
    
    cp2 =  1.0e+03 *[ 0.7333   -2.2127    0.7484];
    cp1 = 1e3*[0.1651   -2.3655   -0.3258];
    set(gca,'cameraposition',cp1 - (cp1-cp2)*timetT/100 );
    %camorbit(-8,-162); %37.5000-(timetT)*37.5/100,-30+(timetT)*30/100)
    surfhandles=findobj('type','surf');
    get(surfhandles,'ambientStrength');
    camlight('head');
    lighting phong
    drawnow;
    h = getframe(G);
    imwrite(imresize(h.cdata, [660 860]),['prefraime',num2str(timetT),'.png']);
% pause(0.5)
        close gcf
end


%% during
for timeT = 1:5:size(smoothed_timetracesX_cooked,2)
    G = figure('Position',[110 110 860 660],'Color',[1 1 1]);
    for i = 1:numel(coo_x)
        surf(x*0.1+400, y*0.1+400, z*0.1+400,ones(size(z))*600); hold on;
        if ~isnan(coo_x(i))
            activity = min(max(smoothed_timetracesX_cooked(i,timeT)*100,0),600);
%             colorX = [0.66 0.66 0.66]*(1-activity/300) + [1 0 0]*(activity/300);
            r = 2;
            surf(x*r+coo_x(i), y*r+coo_y(i), z*r+coo_z(i),ones(size(z))*activity/1); hold on;
        end
    end
    colormap(cmap);
    caxis([-100 600])
    alpha color
    alphamap('increase')
    colorbar;
    shading(gca,'interp')
    hold off;
    axis([min(coo_x) max(coo_x) min(coo_y) max(coo_y) min(coo_z) max(coo_z)])
%     cp2 = cp2;
    set(gca,'cameraposition',cp2);
    surfhandles=findobj('type','surf');
    get(surfhandles,'ambientStrength');
    camlight('head');
    lighting phong
    drawnow;
    h = getframe(G);
    imwrite(imresize(h.cdata, [660 860]),['frame',num2str(timeT),'.png']);
    close gcf
end
for timetT = 1:20
    G = figure('Position',[ 110 110 860 660],'Color',[1 1 1]);
    surf(x*0.1+400, y*0.1+400, z*0.1+400,ones(size(z))*600); hold on;
    surf(x*0.1+400, y*0.1+400, z*0.1+400,ones(size(z))*-100); hold on;
    for i = 1:numel(coo_x)
    if ~isnan(coo_x(i))
    r = 2;
    surf(x*r+coo_x(i), y*r+coo_y(i), z*r+coo_z(i),ones(size(z))*(600-timetT*40)); hold on;
    end
    end
    colormap(cmap);
    caxis([-100 600])
    alpha color
    colorbar;
    shading(gca,'interp')
    hold off;
    axis([min(coo_x) max(coo_x) min(coo_y) max(coo_y) min(coo_z) max(coo_z)])
    %     pause(0.1);
    set(gca,'cameraposition',cp2);
    surfhandles=findobj('type','surf');
    get(surfhandles,'ambientStrength');
    camlight('head');
    lighting phong
    drawnow;
    h = getframe(G);
    imwrite(imresize(h.cdata, [660 860]),['postflame',num2str(timetT),'.png']);
    close gcf
end


% cp = get(gca,'cameraposition');

