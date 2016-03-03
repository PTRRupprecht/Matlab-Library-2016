
MatFileList = dir('Extracted*.mat');
for yy = 1:numel(MatFileList)
    load(MatFileList(yy).name);
    clear fluoXT fluoX
    for trailx = 1:8
        fluoX = [];
        for jj = 1:4
            fluotrace = plane{jj}.timetraces_raw;
            indizes = find(~isnan(sum(fluotrace(:,:,1),1)));
            fluotrace = fluotrace(:,indizes,:);
            for kk = 1:size(fluotrace,2)
                ftrace = smooth(fluotrace(:,kk,trailx),25); F0 = min(ftrace(25:end-25));
                fluotraceX = (fluotrace(:,kk,trailx) - F0)/F0*100;
                fluoX = [fluoX, fluotraceX];
            end
        end
        fluoXT(:,:,trailx) = fluoX;
    end
    figure(3838);
    subplot(1,3,yy); imagesc(corr(squeeze(mean(fluoXT(104:144,:,:),1))),[0 1])
    figure(1818); 
    for trailx = 1:8
        subplot(3,8,trailx + (yy-1)*8); imagesc((1:351)/7.5,[],mean(fluoXT(:,:,trailx),3)',[0 200]);
    end
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
