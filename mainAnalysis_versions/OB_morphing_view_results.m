

MatFileList = dir('Extracted*.mat');
load(MatFileList(1).name);
clear fluoXT fluoX
for trailx = 1:numel(plane{1}.timetraces)
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

subplot(1,numel(MatFileList),yy); imagesc(corr(squeeze(mean(fluoXT(104:148,1:end,:),1))),[0 1])
%     for trailx = 1:size(fluoXT,3)
%         subplot(3,size(fluoXT,3),trailx + (yy-1)*8); imagesc((1:351)/7.5,[],mean(fluoXT(:,:,trailx),3)',[0 200]);
%     end
figure(1218+yy);
%     figure(1299); 
for trailx = 1:size(fluoXT,2)
    subplot(ceil(size(fluoXT,2)/floor(sqrt(size(fluoXT,2)))),floor(sqrt(size(fluoXT,2))),trailx); imagesc((1:351)/7.5,[],squeeze(fluoXT(:,trailx,:))',[0 200]);
    title(strcat('Plane',32,num2str(planesIX(trailx)),', neuron',32,num2str(planesIXX(trailx))));
end

figure(12319); plane_nb = 4; neuron_nb = 6;
windowsize = 30;
ROIX = squeeze(plane{plane_nb}.ROI_map(1,:,:));
[x,y] = find(ROIX == neuron_nb); x = round(mean(x)); y = round(mean(y));
xxx = max(1,x-windowsize):min(512,x+windowsize); yyy = max(1,y-windowsize):min(477,y+windowsize);
ROIXX = ROIX(xxx,yyy);
subplot(4,4,1); imagesc(ROIXX==neuron_nb);
for zz = 1:size(plane{plane_nb}.timetraces_raw,3); subplot(4,4,1+zz); imagesc(plane{plane_nb}.anatomy(xxx,yyy,zz),[-30 70]); end; colormap(gray)
for zz = 1:size(plane{plane_nb}.timetraces_raw,3); subplot(4,4,1+zz); imagesc(plane{plane_nb}.DF_reponse(xxx,yyy,zz),[-0.5 2]); end; colormap(jet)
subplot(3,5,1); imagesc(ROIXX==neuron_nb);
for zz = 1:size(plane{plane_nb}.timetraces_raw,3); subplot(3,5,1+zz); imagesc(plane{plane_nb}.anatomy(xxx,yyy,zz),[-30 70]); end; colormap(gray)
%     for zz = 1:size(plane{plane_nb}.timetraces_raw,3); subplot(3,4,1+zz); imagesc(plane{plane_nb}.DF_reponse(xxx,yyy,zz),[-0.5 2]); end; colormap(jet)






% discard neurons

% p_ix = [1];
% n_ix = [11];
% 
% t_ix = [9 10 11 12];
% 
% planeX = plane;
% for j = 1:numel(p_ix)
%     j
%     planeX{p_ix(j)}.ROI_map(planeX{p_ix(j)}.ROI_map == n_ix(j)) = 0;
%     planeX{p_ix(j)}.timetraces_raw(:,n_ix(j),:) = NaN;
%     planeX{p_ix(j)}.timetraces(:,n_ix(j),:) = NaN;
%     planeX{p_ix(j)}.timetraces(:,n_ix(j),:) = NaN;
% end
% for k = 1:4
%     planeX{k}.timetraces(:,:,t_ix) = [];
%     max_timetrace = size(planeX{k}.timetraces_raw(:,:,t_ix),2);
%     planeX{k}.timetraces(:,(max_timetrace+1):end,:) = [];
%     planeX{k}.timetraces_raw(:,:,t_ix) = [];
%     planeX{k}.meta.numberframes(t_ix) = [];
%     planeX{k}.ROI_map(t_ix,:,:) = [];
%     planeX{k}.DF_reponse(:,:,t_ix) = [];
%     planeX{k}.anatomy(:,:,t_ix) = [];
% end
% plane_04_03_His = planeX;
% 
% save('Fish2_15_03_Trp.mat','plane')



