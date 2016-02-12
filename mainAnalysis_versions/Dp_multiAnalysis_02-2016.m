clear all
global clut2b timetracesX_X ROI_map_X movie_AVG_X
load clut2b
% load list of files
% cd distorted
FileList_0 = dir('Test_028*_.tif');
% cd ..
FileList = FileList_0;%dir('OB_4planes_paradigm1*_undistorted.tif');
FileListX = FileList;
for k = 1:numel(FileList)
    FileList(k) = FileListX(numel(FileList)+1-k);
end
% cd distorted
% read metadata
clear meta
[A,result,meta.framerate,meta.zstep,meta.zoom,meta.motorpositions,meta.scalingfactors] = read_metadata_function(FileList_0(1).name);
% cd ..
for kkk = 1:numel(FileList)
    L{kkk} = imfinfo(FileList(kkk).name);
    meta.height = L{kkk}(1).Height;
    meta.width = L{kkk}(1).Width;
    meta.numberframes(kkk) = numel(L{kkk});
end
% initialize 3D matrix
binning = 1;

meta.framerate = meta.framerate/binning;
nb_frames_total = sum(floor(meta.numberframes/binning));
clear movie movie_p
useful_range_start = 60;

nb_planes = 5;
meta.framerate = meta.framerate/nb_planes;
nb_frames_perplane = floor(nb_frames_total/nb_planes);
for pp = 1:nb_planes
    
    
    
    counter_planes{pp} = 0;
    movie_p{pp} = zeros(meta.height,meta.width,nb_frames_perplane);
    for kkk = 1:numel(FileList)
        kkk
        % read raw data
        indicator = rem(sum(meta.numberframes(1:kkk-1)),nb_planes);
        if indicator < pp
            startingpoint = pp - indicator;
        else
            startingpoint = nb_planes + pp - indicator;
        end
        nb_frames_this_time = ceil((meta.numberframes(kkk)-startingpoint+1)/nb_planes);
        [movie_p{pp}(:,:,((counter_planes{pp}+1):(counter_planes{pp}+nb_frames_this_time))),movie_AVG_X{kkk,pp}] = read_movie(FileList(kkk).name,meta.width,meta.height,nb_frames_this_time,startingpoint,binning,L{kkk},nb_planes);
        counter_planes{pp} = counter_planes{pp} + nb_frames_this_time;
    end
    
    temp = movie_p{pp};
    movie_p{pp}(:,:,1:1600) = temp(:,:,801:end);
    movie_p{pp}(:,:,1601:2400) = temp(:,:,1:800);
    
    % discard einschwingvorgang
    movie_p{pp} = movie_p{pp}(:,:,useful_range_start:end);
    % align movie intra-trial
    reference{pp} = mean(movie_p{pp}(:,:,(-200:200) + round(nb_frames_total/2/nb_planes)),3);
    [movie_p{pp},offsety_resolved,offsetx_resolved] = alginWithinTrial_revised(movie_p{pp},reference{pp});
    movie_AVG_X{pp} = mean(movie_p{pp},3);
    % check alignment visually
    % implay(movie_p{1}(:,:,1:20:end)/max(movie_p{1}(:)))
    % estimate offset of PMTs
    offset = quantile(movie_AVG_X{1}(:),0.02);
    % dF over F
end

% correct for DAQ board transient
% template = zeros(size(movie_p{pp}(:,:,1)));
% for pp = 1:nb_planes
%     dummy = mean(movie_p{pp}(:,:,1:end),3);
%     template(1:2:end,:) = template(1:2:end,:) + repmat((  nanmean(dummy(1:2:end,:),1)  ) , [size(template,2)/2 1]);
%     template(2:2:end,:) = template(2:2:end,:) +  repmat(( nanmean(dummy(2:2:end,:),1)), [size(template,2)/2 1]);
% %     imagesc(template)
% end
% template = template/nb_planes;

% for pp = 1:nb_planes
%     movie_p{pp} = movie_p{pp} + repmat(template,[1 1 size(movie_p{pp},3)]);
% end

% LL = zeros(size(movie_AVG_X{pp}));
for pp = 1:nb_planes
    movie_AVG_X{pp} = mean(movie_p{pp},3);
%     movie_AVG_X{pp} = movie_AVG_X{pp} - template;
% LL =  LL+ movie_AVG_X{pp};
end

for pp = 1:nb_planes
    plot1 = 34; plot2 = 0; DF_movie_yesno = 0; % figure number
    [DF_reponse{pp},DF_master{pp},DF_movie] = dFoverF(movie_p{pp},offset,meta.framerate,plot1,plot2,DF_movie_yesno);
    drawnow;
    % local correlation map (computational slightly expensive, but good
    % alternative to dF over F map; movie has to be 2^x for width and height
    tilesize = 16;
    localCorrelations{pp} = localCorrelationMap(movie_p{pp},tilesize);
    % localCorrelations = zeros(size(DF_reponse));
end

figure(88);
colormap(paruly)
for pp = 1:5
    anatomy{pp} = undistort_stack(mean(movie_p{pp}(:,:,800:1300),3),meta.zoom);
    subplot(2,3,pp); imagesc(adapthisteq((anatomy{pp}-min(anatomy{pp}(:)))/(max(anatomy{pp}(:))-min(anatomy{pp}(:))),'NumTiles',[8 8])); colormap(gray)
end



% FOV changes
xd = [-98.5:29.5:150];
[1-0.1/100*xd]
1-0.3/100*xd



% figure(922);
% for pp = 1:4
%     F0_window = round([1 23]*meta.framerate); F0_window = F0_window(1):F0_window(2);
%     response_window1 = round([27.7 58]*meta.framerate); response_window1 = response_window1(1):response_window1(2);
%     response_window2 = round([97.7 128]*meta.framerate); response_window2 = response_window2(1):response_window2(2);
%     filtersize = 2;
%     kernelX = fspecial('gaussian',[filtersize filtersize], 2);
%     F0 = mean(movie_p{pp}(:,:,F0_window),3)+80;
%     F0 = conv2(F0,kernelX,'same');
%     DF1 = mean(movie_p{pp}(:,:,response_window1),3);
%     DF1 = (DF1 - F0)./(F0-offset);
%     DF1_reponse = conv2(DF1,kernelX,'same');
%     DF2 = mean(movie_p{pp}(:,:,response_window2),3);
%     DF2 = (DF2 - F0)./(F0-offset);
%     DF2_reponse = conv2(DF2,kernelX,'same');
%     %     DF1_reponse = undistort_stack(DF1_reponse,meta.zoom,0);
%     %     DF2_reponse = undistort_stack(DF2_reponse,meta.zoom,0);
%     DF1X{pp} = DF1_reponse;
%     DF2X{pp} = DF2_reponse;
%     subplot(2,4,pp); imagesc(DF1_reponse,[-0.660 1.9]); axis off equal
%     subplot(2,4,pp+4); imagesc(DF2_reponse,[-0.66 1.9]); axis off equal
% end


for pp = 1:5

    
   pp = 5
    df_scale = [-10 200];
    last_ii = 1;
    % extract cellular time traces . semi-automated ROI-detection
    ROI_map_input = ROI_mapXX{pp};
    [ROI_mapX,timetracesX,timetracesX_raw] = timetraces_singleplane(movie_p{pp},movie_AVG_X{pp},offset,DF_reponse{pp},DF_master{pp},localCorrelations{pp},df_scale,ROI_map_input,meta,last_ii,movie_AVG_X{pp});
    % figure, imagesc(conv2(timetracesX,fspecial('gaussian',[25 1],23),'same'));
    ROI_map_input = ROI_mapX;


    ROI_mapXX{pp} = ROI_mapX;
    timetracesX_X{pp} = timetracesX;
    timetracesX_X_raw{pp} = timetracesX_raw;

    [size(timetracesX_X_raw{1},2), size(timetracesX_X_raw{2},2), size(timetracesX_X_raw{3},2),size(timetracesX_X_raw{4},2),size(timetracesX_X_raw{5},2)]
    sum([size(timetracesX_X_raw{1},2), size(timetracesX_X_raw{2},2), size(timetracesX_X_raw{3},2),size(timetracesX_X_raw{4},2),size(timetracesX_X_raw{5},2)])
end
    
%% plotting anatomy
figure(811);
colormap(lines)
for pp = 1:5
    ROI_mapXXZ{pp} = undistort_stack_integers(ROI_mapXX{pp},meta.zoom);
    subplot(5,1,pp); imagesc(ROI_mapXXZ{pp}); axis off equal; colormap(lines)
end


%% concenating all timetraces
fullX =     [(timetracesX_X_raw{1}), (timetracesX_X_raw{2}), (timetracesX_X_raw{3}),(timetracesX_X_raw{4}),(timetracesX_X_raw{5})];

fullX = [];
for pp = 1:5
    fullX = [fullX, timetracesX_X_raw{pp}(:,find(ppClusterX{pp} == 2))];
end

[COEFF,SCORE] = princomp(fullX); % PCA, for no special reason

figure(77),
for i = 1:3
    subplot(3,1,i); plot(SCORE(:,i));
end
%% select useful timetraces

global clusterX
for pp = 1:5
    
    for k = 1:size(timetracesX_X_raw{pp},2)
        handle1 = figure(2198); plot(smooth(timetracesX_X_raw{pp}(:,k),5));
        set(gcf, 'units','normalized','Position',[0.1 0.1 0.55 0.620]);
        set(gcf, 'WindowKeyPressFcn', {@chooseCluster,k,handle1});
        grid on;
        waitfor(gcf);
    end
    ppClusterX{pp} = clusterX;
end
save('ClusterX.mat','ppClusterX')


%% coloring ROIs according to their correlation with cluster centers in RGB
cmap = jet(999);
cmap(500,:) = [ 1 1 1];
figure(811);
for pp = 1:5
    pp
    D = repmat(ROI_mapXX{pp},[1 1 3]);
    for j = 1:max(D(:))
        [x,y] = find(ROI_mapXX{pp} == j);
        if ~isempty(x)
            A1 = corr(timetracesX_X_raw{pp}(:,j),G(:,2));
            A2 = corr(timetracesX_X_raw{pp}(:,j),G(:,3));
            A3 = corr(timetracesX_X_raw{pp}(:,j),G(:,4));
            for k = 1:numel(x)
                D(x(k),y(k),:) = cat(3,((A1)),((A2)),((A3)));
            end
        end
    end
%     D(D == 0) = 1;
    D(D == 0) = NaN;
    D = D + 0.4559;
%     P(pp) = min(D(:));
%     R(pp) = max(D(:));
    D = D/(0.8356+0.4559);
    
    D(isnan(D(:))) = [0.5];
    D = undistort_stack_integers(D,meta.zoom);
    subplot(2,3,pp); imagesc(D,[0 1]); axis off equal; %colormap(cmap)
end

figure(441); hold on; 
plot(G(:,2),'r');
plot(G(:,3)+300,'g');
plot(G(:,4)+600,'b');



%% correlation maps in RGB !!!!

for k = 1:3
    timetrace_X = G(:,k+1);
    for pp = 1:5
        k
        pp
        xcorr_tot{k,j} = zeros(512);
        for j = 1:512
            if mod(j,50) == 0; disp(j/size(movie_p{pp},1)); end
            xcorr_tot{k,pp}(j,:) = corr(timetrace_X,squeeze(movie_p{pp}(j,:,:))');
        end
    end
end
cmap = jet(999);
cmap(500,:) = [ 1 1 1];
figure(8311);
for pp = 1:5
    D = cat(3,xcorr_tot{1,pp},xcorr_tot{2,pp},xcorr_tot{3,pp});
    D = D + 0.3986;
%     P(pp) = min(D(:));
%     R(pp) = max(D(:));
    D = D/(0.6809+0.3985);
    D = D*1.5; D = D - 0.25;
    D = undistort_stack(D,meta.zoom);
    subplot(2,3,pp); imagesc(max(min(D,1),0)); axis off equal; %colormap(cmap)
end


%% plotting the yz-scanning trajectory

tz = linspace(98.5+15,-150-15,256*9);
% tz = tz(end:-1:1);

ty = (-128:127)*0.7258;
ty = [ty,ty,ty,ty,ty,ty,ty,ty,ty].*(1-0.1./100.*tz);


figure, hold on;
for k = 1:9
    plot(ty(((k-1)*256+1):k*256),tz(((k-1)*256+1):k*256));
end
xd = [-98.5:29.5:150];
1-0.1/100*xd
1-0.3/100*xd




%% sorting and clustering ...


FF = fullX;

nb_clusters = 4; % choose number of clusters to detect

detect_traces = find(nansum(FF));
dF = FF(:,detect_traces);
corrMatrix = corr(dF);
corrMatrixOld = corrMatrix;
corrMatrix = corrMatrix - eye(size(corrMatrix));
dissimilarity = 1 - corrMatrix';

% perform complete linkage clustering
Z = linkage(dissimilarity,'complete');
% decide on a cutoff
cutoff = 10;
l_groups = 1;
while l_groups < nb_clusters
    cutoff = cutoff - 0.1;
    groups = cluster(Z,'cutoff',cutoff,'criterion','distance');
    l_groups = length(unique(groups));
end

[XI,IX] = sort(groups);
dF_ordered = dF(:,IX);

figure(7); hold on;
for k = 1:nb_clusters
    indizet = find(XI == k);
    plot(mean(dF_ordered(:,indizet),2)+300*k)
    G(:,k) = mean(dF_ordered(:,indizet),2);
end
hold off;
corrMatrix2 = corr(dF_ordered);

% figure, imagesc(dF_ordered)
figure(8), imagesc(corrMatrix2)
% [IX,XI] = sort(max(dF_ordered));

cmap = distinguishable_colors(size(dF_ordered,2));
figure(4); hold on;
for j = 1:size(dF_ordered,2)
    plot(smooth(dF_ordered(:,IX(end-j-1)) + 440*j,4),'Color',cmap(j,:));
end
    

%% plotting some timetraces with various coloring
% fullX = [];
figure(111), hold on;
counter = 0;
cmap = distinguishable_colors(5);
for pp = 1:5
%     fullX = [fullX, timetracesX_X_raw{pp}(:,find(ppClusterX{pp} == 2))];
    F = timetracesX_X_raw{pp}(:,find(ppClusterX{pp} == 2));
    for k = 1:size(F,2)
        if sum(F(:,k)) < 0; F(:,k) = F(:,k)*-1; end
        F0 = min(smooth(F(:,k),30));
        if F0 < 40; F0 = 40; end
        F(:,k) = (F(:,k)-F0)/(F0)*100;
        if sum(F(:,k)) < 0; F(:,k) = F(:,k)*-1; end
        
%         F(:,k) = F(:,k) - min(smooth(F(:,k),20));
    end
    
    [~, strongest] = sort(max(F),'descend');
    indises = strongest(1:min(15,numel(strongest)));
    for j = 1:numel(indises)
        tt = F(:,indises(numel(indises)+1-j));
        A1 = corr(tt,G(:,2));
        A2 = corr(tt,G(:,3));
        A3 = corr(tt,G(:,4));
        
        plot((1:2341)/6,smooth(F(:,indises(numel(indises)+1-j)) + 640*counter,2),'Color',([0 0 0 ])); %[A1 A2 A3]+ 0.4559)/(0.8356+0.4559)); %cmap(pp,:));
        counter = counter + 1;
    end
end


%% okay ... not so useful any more!

DP_total.ROI_mapXX = ROI_mapXX;
DP_total.timetracesX_X = timetracesX_X;
DP_total.timetracesX_X_raw = timetracesX_X_raw;
DP_total.anatomy = anatomy;
DP_total.DF_reponse = DF_reponse;
DP_total.DF_master = DF_master;
DP_total.localCorrelations = localCorrelations;
save(strcat(FileList(1).name(1:end-4),'_timetrace.mat'),'DP_total','meta','offset');

