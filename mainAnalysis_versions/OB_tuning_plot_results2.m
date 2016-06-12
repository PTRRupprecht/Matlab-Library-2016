


MatFileList = dir('Extracted_Data_*.mat');
for yy = 1:numel(MatFileList)
    load(MatFileList(yy).name);
    clear fluoXT fluoX
    fluoXT = [];
    for trailx = 1:numel(plane{1}.timetraces)
        fluoX = [];
        planesIX = [];
        planesIXX = [];
        for jj = 1:numel(plane)
            fluotrace = plane{jj}.timetraces_raw{trailx};
%             indizes = find(~isnan(sum(fluotrace(:,:,1),1)));
%             fluotrace = fluotrace(:,indizes);
            fluotrace = fluotrace(:,:);
            for kk = 1:size(fluotrace,2)
                ftrace = smooth(fluotrace(:,kk),25); F0 = min(ftrace(25:end-25));
                fluotraceX = (fluotrace(:,kk) - F0)/F0*100;
                fluoX = [fluoX, smooth(fluotraceX,5)];
                planesIX = [planesIX,jj];
%                 planesIXX = [planesIXX,indizes(kk)];
                planesIXX = [planesIXX,kk];
            end
        end
        fluoXT{trailx} = fluoX;
    end
        
    figure(760);
    for trailx = 1:numel(plane{1}.timetraces)
        subplot(2,3,trailx); imagesc(fluoXT{trailx}')
    end
    
    
    
    [corrMap,timepoints] = corrNaNs(fluoXT);
    
    select_points = [1 2 13 14 15 16];
    select_points = [1 2 4 6 7 9 10 12 13:16];
    select_points = [1:6];
    select_times = [];
    for j = 1:numel(select_points); select_times = [select_times, (1+sum(timepoints(1:select_points(j)))-timepoints(select_points(j))):sum(timepoints(1:select_points(j)))]; end
    
    figure(9891), hold on; for k = 1:6; plot(nanmean(fluoXT{k},2)); end
    figure(127+yy);
    subplot(1,1,1); imagesc(corrMap(select_times,select_times)); axis equal; %axis([0  size(corrMap,1) 0 size(corrMap,1)])
end
    
    for kk = 1:numel(select_points)
        for ii = 1:numel(select_points)
            timek = sum(timepoints(1:(select_points(kk)-1))) + 1;
            timei = sum(timepoints(1:(select_points(ii)-1))) + 1;
            timekk = sum(timepoints(1:(select_points(kk)-1))) + timepoints(select_points(kk));
            timeii = sum(timepoints(1:(select_points(ii)-1))) + timepoints(select_points(ii));
            for jj = 1:min(numel(timei:timeii),numel(timek:timekk))
                Corr3D(kk,ii,jj) = corrMap( timek+jj-1,timei+jj-1);
            end
        end
    end
    
    FF = reshape(permute(Corr3D,[3 1 2]),[551 144]); 
    FX = FF;
    FX(:,[1:13:end]) = [];
    figure, imagesc((1:551)/7.5,[],FX')
    figure, plot(mean(FX,2))
end