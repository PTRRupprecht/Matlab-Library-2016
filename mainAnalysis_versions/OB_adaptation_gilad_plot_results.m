
MatFileList = dir('Fish*X.mat');
for yy = 1%numel(MatFileList)
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
    figure(3238);
    subplot(1,1,yy); imagesc(corr(squeeze(mean(fluoXT(104:148,1:end,:),1))),[0 1])
%     for trailx = 1:size(fluoXT,3)
%         subplot(3,size(fluoXT,3),trailx + (yy-1)*8); imagesc((1:351)/7.5,[],mean(fluoXT(:,:,trailx),3)',[0 200]);
%     end
    figure(1218+yy);
%     figure(1299); 
    for trailx = 1:size(fluoXT,2)
        subplot(ceil(size(fluoXT,2)/floor(sqrt(size(fluoXT,2)))),floor(sqrt(size(fluoXT,2))),trailx); imagesc((1:351)/7.5,[],squeeze(fluoXT(:,trailx,:))',[0 200]);
        title(strcat('Plane',32,num2str(planesIX(trailx)),', neuron',32,num2str(planesIXX(trailx))));
    end
    
    figure(12319); plane_nb = 3; neuron_nb = 9;
    windowsize = 20;
    ROIX = squeeze(plane{plane_nb}.ROI_map(1,:,:));
    [x,y] = find(ROIX == neuron_nb); x = round(mean(x)); y = round(mean(y));
    xxx = max(1,x-windowsize):min(512,x+windowsize); yyy = max(1,y-windowsize):min(477,y+windowsize);
    ROIXX = ROIX(xxx,yyy);
    subplot(4,4,1); imagesc(ROIXX==neuron_nb);
    for zz = 1:size(plane{plane_nb}.timetraces_raw,3); subplot(4,4,1+zz); imagesc(plane{plane_nb}.anatomy(xxx,yyy,zz),[-30 70]); end; colormap(gray)
    for zz = 1:size(plane{plane_nb}.timetraces_raw,3); subplot(4,4,1+zz); imagesc(plane{plane_nb}.DF_reponse(xxx,yyy,zz),[-0.5 2]); end; colormap(jet)

end











% discard neurons

p_ix = [];
n_ix = [];

t_ix = [9 10];

planeX = plane;
for j = 1:numel(p_ix)
    j
    planeX{p_ix(j)}.ROI_map(planeX{p_ix(j)}.ROI_map == n_ix(j)) = 0;
    planeX{p_ix(j)}.timetraces_raw(:,n_ix(j),:) = NaN;
    planeX{p_ix(j)}.timetraces(:,n_ix(j),:) = NaN;
end
for k = 1:4
    planeX{k}.timetraces(:,:,t_ix) = [];
    planeX{k}.timetraces_raw(:,:,t_ix) = [];
    planeX{k}.meta.numberframes(t_ix) = [];
    planeX{k}.ROI_map(t_ix,:,:) = [];
    planeX{k}.DF_reponse(:,:,t_ix) = [];
    planeX{k}.anatomy(:,:,t_ix) = [];
end
plane_04_03_His = planeX;

save('Fish5_03_03_His.mat','plane_04_03_His')




FileLischt = dir('Fish*.mat');
for j = 1:numel(FileLischt)
    j = j + 1
    F = load(FileLischt(j).name);
    keyboard
    if 0
        plane = F.plane_04_03_His
    end
    for i = 1:4
        plane{i}.timetraces(:,(size(plane{i}.timetraces_raw,2)+1):end,:) = [];
    end
%     timecourse = [];
%     for i = 1:4
%     timecourse = [timecourse,(nanmean(plane{i}.timetraces,3))];
%     end
%     timecourse = nanmean(nanstd(plane{i}.timetraces,0,3),2);
%     figure(88); plot((timecourse))
%     figure(88); imagesc(squeeze(timecourse))
    save(strcat(FileLischt(j).name(1:end-4),'X.mat'),'plane');
end


    
    
        if size(plane{i}.timetraces,2)
    
    timecourse = 



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
