
function [firingZ,EM] = timingMap(timetraces,ROI_map_now,framerate,valve_switch_sec,smoothing1,smoothing2,showMap)

    A = timetraces;
    
    nbtiles = ceil(size(timetraces,2)/100);
    clear firingXX
    for ii = 1:nbtiles
        AX = A(:,((ii-1)*100+1):min((ii-1)*100+100,max(ROI_map_now(:)))); %(framerate*round(valve_switch_sec):framerate*round(valve_switch_sec+20),:);
        clear SX firingX
        for i = 1:size(AX,2)
            B = AX(:,i);

            S = B - circshift(B,smoothing1);  S = S((smoothing1+1):(end-smoothing1));
            Kuckuck = smooth(S,smoothing2); % derivative, smoothed
            positions = strfind(Kuckuck'>std(Kuckuck)*1.0,ones(smoothing1-2,1)'); % find steadily rising slope
            if isempty(positions)
                firing = 0;
            else
                firing = positions(1) + 0;
                La = smooth(B,5); wiring = La(firing + smoothing1);
                if wiring < 25
                    firing = 0;
                end
            end;
            firingX(i) = firing + smoothing1;
            SX(:,i) = Kuckuck;
        end

        % figure(53), imagesc(SX,[0 100]); hold on; plot(firingX,'y.','MarkerSize',32); hold off
        % figure(6), imagesc(AX,[0 100]); hold on; plot(firingX+smoothing1,'y.','MarkerSize',32); hold off

        endofdata = find(AX(1,:)==0,1,'first') - 1;
        if ~isempty(endofdata); AX = AX(:,1:endofdata); end;
        AX =  conv2( AX, fspecial('gaussian',[10 1],6),'same');
        figure(53), imagesc(AX,[0 100]); hold on; plot(firingX,'y.','MarkerSize',32); hold off;
        set(gcf, 'units','normalized','Position',[1.1 0.1 0.8 0.8]); 
        figure(54); set(gcf, 'units','normalized','Position',[0.02 0.55 0.25 0.35]); 
        while ishandle(54)
            figure(54); try; plot(AX(:,coordinates(1))); hold on; plot( coordinates(2),A( coordinates(2),coordinates(1)),'.k','MarkerSize',32); hold off; end
            figure(53); imagesc(AX,[0 100]); hold on; plot(firingX,'y.','MarkerSize',32); hold off; set(gcf, 'units','normalized','Position',[0.3 0.3 0.6 0.6]); 
            
            try
                l = ginput(1);
            catch
                close 54
            end
            coordinates = round(l);
            firingX(coordinates(1)) = coordinates(2);
        end
    
        firingXX( ((ii-1)*100+1):min((ii-1)*100+100,max(ROI_map_now(:))) ) = firingX;
        
    end

    firingZ = (firingXX-1)/framerate;

    All = ROI_map_now;
    EM = zeros(size(All));
    for i = 1:max(All(:))
        EM(All == i) = firingZ(i);
    end
    if showMap
        figure(3322), imagesc(EM,[2 5]); colormap(jet)
    end

end
