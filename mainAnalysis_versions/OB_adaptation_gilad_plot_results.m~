
MatFileList = dir('Extracted*.mat');
for yy = 1:numel(MatFileList)
    load(MatFileList(yy).name);
    clear fluoXT fluoX
    for trailx = 1:size(plane{1}.timetraces,3)
        fluoX = [];
        planesIX = [];
        planesIXX = [];
        for jj = 1:numel(plane)
            fluotrace = plane{jj}.timetraces_raw;
            indizes = find(~isnan(sum(fluotrace(:,:,1),1)));
            fluotrace = fluotrace(:,indizes,:);
            for kk = 1:size(fluotrace,2)
                ftrace = smooth(fluotrace(:,kk,trailx),25); F0 = min(ftrace(25:end-25));
                fluotraceX = (fluotrace(:,kk,trailx) - F0)/F0*100;
                fluoX = [fluoX, smooth(fluotraceX,5)];
                planesIX = [planesIX,jj];
                planesIXX = [planesIXX,indizes(kk)];
            end
        end
        fluoXT(:,:,trailx) = fluoX;
    end
    figure(3838);
    subplot(1,numel(MatFileList),yy); imagesc(corr(squeeze(mean(fluoXT(104:148,1:end,:),1))),[0 1])
%     for trailx = 1:size(fluoXT,3)
%         subplot(3,size(fluoXT,3),trailx + (yy-1)*8); imagesc((1:351)/7.5,[],mean(fluoXT(:,:,trailx),3)',[0 200]);
%     end
    figure(1818+yy); 
    for trailx = 1:size(fluoXT,2)
        subplot(ceil(size(fluoXT,2)/floor(sqrt(size(fluoXT,2)))),floor(sqrt(size(fluoXT,2))),trailx); imagesc((1:351)/7.5,[],squeeze(fluoXT(:,trailx,:))',[0 200]);
        title(strcat('Plane',32,num2str(planesIX(trailx)),', neuron',32,num2str(planesIXX(trailx))));
    end
    
%     figure(19319); plane_nb = 1; neuron_nb = 10;
%     ROIX = squeeze(plane{plane_nb}.ROI_map(1,:,:));
%     [x,y] = find(ROIX == neuron_nb); x = round(mean(x)); y = round(mean(y));
%     xxx = max(1,x-40):min(512,x+40); yyy = max(1,y-40):min(477,y+40);
%     ROIXX = ROIX(xxx,yyy);
%     subplot(3,4,1); imagesc(ROIXX==neuron_nb);
%     for zz = 1:10; subplot(3,4,1+zz); imagesc(plane{plane_nb}.anatomy(xxx,yyy,zz),[-30 70]); end; colormap(gray)

end



% figure(24);
% for jj = 1:4
%     fluotrace = plane{jj}.timetraces_raw;
%     indizes = find(~isnan(sum(fluotrace(:,:,1),1)));
%     fluotrace = fluotrace(:,indizes,:);
%     for kk = 1:size(fluotrace,2)
%         for xx = 1:size(fluotrace,3)
%             ftrace = smooth(fluotrace(:,kk,xx),25); F0 = min(fluotrace(25:end-25));
%             fluotrace(:,kk,xx) = (fluotrace(:,kk,xx) - F0)/F0*100;
%         end
%     end
%     for j = 1:8
%         subplot(4,8,j+(jj-1)*8); imagesc([],(1:351)/7.5,fluotrace(:,:,j),[0 600]);
%     end
% end
