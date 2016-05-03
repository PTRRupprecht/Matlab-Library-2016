
%% load data (plane, morphing paradigms, Filenames)
MatFileList = dir('Extracted*fish1_*.mat');
load(MatFileList(1).name);
load('paradigms.mat');
F = dir('Fish1_*.tif');
FileNames = cat(1,F.name);


%% load paradigms
pdgLUT{1} = paradigm.p91; pdgLUT{2} = paradigm.p92;
pdgLUT{3} = paradigm.p93; pdgLUT{4} = paradigm.p94;
pdgLUT{5} = paradigm.p95; pdgLUT{6} = paradigm.p96;
try; pdgLUT{7} = paradigm.p97; end
try; pdgLUT{8} = paradigm.p98; end
nb_trials = size(plane{1}.anatomy,3);
smoothing = 5;

%% make anatomy clickable (RBG overlay image)
all_anatomy = zeros(1024);
select_trial = 7;
all_anatomy(1:512,1:512) = plane{1}.anatomy(:,:,select_trial); all_ROI(1:512,1:512) = squeeze(plane{1}.ROI_map(select_trial,:,:));
all_anatomy(513:1024,1:512) = plane{2}.anatomy(:,:,select_trial); all_ROI(513:1024,1:512) = squeeze(plane{2}.ROI_map(select_trial,:,:));
all_anatomy(1:512,513:1024) = plane{3}.anatomy(:,:,select_trial); all_ROI(1:512,513:1024) = squeeze(plane{3}.ROI_map(select_trial,:,:));
all_anatomy(513:1024,513:1024) = plane{4}.anatomy(:,:,select_trial); all_ROI(513:1024,513:1024) = squeeze(plane{4}.ROI_map(select_trial,:,:));
planeColor(1:512,1:512) = 1; 
planeColor(513:1024,1:512) = 2;
planeColor(1:512,513:1024) = 3; 
planeColor(513:1024,513:1024) = 4;
% RGB image
all_anatomy(all_anatomy> quantile(all_anatomy(:),0.995)) = quantile(all_anatomy(:),0.995);
all_anatomy = all_anatomy - min(all_anatomy(:));
all_anatomy = all_anatomy/max(all_anatomy(:));
scaling = 2.5;
ROI_map_thresh = all_ROI; ROI_map_thresh(ROI_map_thresh>1) = 1;
RGBimage = repmat(all_anatomy*scaling,[1 1 3]); 
RGBimage(:,:,1:2) = RGBimage(:,:,1:2) + repmat(ROI_map_thresh/scaling,[1 1 2]);
maxVal = max(RGBimage(:)); RGBimage = RGBimage/maxVal;

%% parameters for callback function
parameters.nb_trials = nb_trials;
parameters.pdgLUT = pdgLUT;
parameters.plane = plane;
parameters.FileNames = FileNames;
parameters.all_ROI = all_ROI;
parameters.planeColor = planeColor;
parameters.smoothing = smoothing;
parameters.paradigms = str2num([FileNames(:,end-9)]); %#ok<ST2NM>
% framerate = 7.5 Hz

%% figure and callback function
figure(918), imagesc(RGBimage); axis equal off
akZoom_PR();
set(gcf, 'WindowKeyPressFcn', {@OB_morphing_selectROI_9x,parameters});


%% discard non-responsive cells
% 2 = super
% 1 = non-responsive
% 4 = intermediate or little informative about stimulus identity

for kk = 1:4
    clear clusterX
    global clusterX
    for i = 1:size(plane{kk}.timetraces_raw{1},2)
        plane_nb = kk;
        neuron_nb = i;
        handle1 = figure(2198);
        nb_y = 4;
        nb_x = 8;
        for k = 1:parameters.nb_trials
            % plot paradigm
            row_pos = floor(k/nb_x-0.001);
            trial_timepoints = size(parameters.plane{plane_nb}.timetraces{k},1);
            trial_timepoints2 = min(floor(trial_timepoints/7.5*100),size(parameters.pdgLUT{parameters.paradigms(k)},1));
            ax(k) = subplot(nb_y,nb_x,k + row_pos*nb_x);
            h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[parameters.pdgLUT{parameters.paradigms(k)}(1:10:trial_timepoints2,2); 0],'r');
            set(h, 'EdgeColor','none', 'FaceAlpha', 0.4); hold on;
            h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[parameters.pdgLUT{parameters.paradigms(k)}(1:10:trial_timepoints2,3); 0],'b');
            set(h, 'EdgeColor','none', 'FaceAlpha', 0.4);
            axis([0 trial_timepoints/7.5 0 3])
            plot((1:trial_timepoints)/7.5+50/7.5,smooth(parameters.plane{plane_nb}.timetraces{k}(1:trial_timepoints,neuron_nb),parameters.smoothing)/70,'k');
            hold off;

            % plot ROIs / neuron anatomy
            windowsize = 30;
            ROIX = squeeze(parameters.plane{plane_nb}.ROI_map(1,:,:));
            [x,y] = find(ROIX == neuron_nb); x = round(mean(x)); y = round(mean(y));
            xxx = max(1,x-windowsize):min(512,x+windowsize); yyy = max(1,y-windowsize):min(477,y+windowsize);
            ax(k+2*parameters.nb_trials) = subplot(nb_y,nb_x,k+ (row_pos)*nb_x+nb_x); imagesc(parameters.plane{plane_nb}.anatomy(xxx,yyy,k),[-30 70]); colormap(gray)
            title(parameters.FileNames(k,:),'Interpreter','None');
        end
        set(gcf, 'units','normalized','Position',[0.1 0.1 0.55 0.620]);
        set(gcf, 'WindowKeyPressFcn', {@chooseCluster,i,handle1});
        waitfor(gcf);
    end
    clusterK{kk} = clusterX;
end

save('clusterK.mat','clusterK')



%% show all responses for each trial, clustered automatically using a reference trial
reference_trial = 7;
tracesX = []; clear TTT;
for trx = 1:size(plane{1}.anatomy,3)
    traces = [];
    traces0 = [];
    numbers = [];
    qualitycheck = [];
    for k = 1:4
        x_length = size(plane{k}.timetraces_raw{reference_trial},2);
        temp = plane{k}.timetraces{trx}(:,1:x_length);
        temp0 = plane{k}.timetraces{reference_trial}(:,1:x_length);
        traces = [traces, temp];
        retro_ix = [k*ones(x_length,1)'; 1:x_length];
        numbers = [numbers, retro_ix];
        qualitycheck = [qualitycheck, clusterK{k}];
        traces0 = [traces0, temp0];
    end
%     idx = find(~isnan(sum(traces0)) & sum(traces0) ~= 0);
    idx = find(qualitycheck == 2);
    nb_clusters = 5;
    [~,XI,IX] = cluster_traces(traces0(:,idx),nb_clusters);
    traces_ordered = traces(:,idx(IX));

    tracesX = [tracesX;traces_ordered];
    
    % only for later investigation for single trials
    TTT{trx} = traces_ordered;
    
    figure(2);
    subplot(4,6,trx), imagesc((1:size(traces_ordered,1))/7.5+50/7.5-1.5,[],traces_ordered',[-40 400]);
    trial_timepoints = size(parameters.plane{1}.timetraces{trx},1);
    trial_timepoints2 = min(floor(trial_timepoints/7.5*100),size(parameters.pdgLUT{parameters.paradigms(trx)},1));
    xlabel('time [sec]');
    subplot(4,6,trx+12),  h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[parameters.pdgLUT{parameters.paradigms(trx)}(1:10:trial_timepoints2,2); 0],'r');
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.4); hold on;
    h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[parameters.pdgLUT{parameters.paradigms(trx)}(1:10:trial_timepoints2,3); 0],'b');
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.4);
    axis([min((1:size(traces_ordered,1))/7.5+50/7.5-1.5) max((1:size(traces_ordered,1))/7.5+50/7.5-1.5) 0 3])

    %     figure(5); 
    %        subplot(2,6,trx),imagesc(corr(traces_ordered));
    % 
    %     figure(4)
    %         trial_timepoints = size(parameters.plane{1}.timetraces{trx},1);
    %         trial_timepoints2 = min(floor(trial_timepoints/7.5*100),size(parameters.pdgLUT{parameters.paradigms(trx)},1));
    %         subplot(2,6,trx),  h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[parameters.pdgLUT{parameters.paradigms(trx)}(1:10:trial_timepoints2,2); 0],'r');
    %         set(h, 'EdgeColor','none', 'FaceAlpha', 0.4); hold on;
    %         h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[parameters.pdgLUT{parameters.paradigms(trx)}(1:10:trial_timepoints2,3); 0],'b');
    %         set(h, 'EdgeColor','none', 'FaceAlpha', 0.4);
end

%% possible switching cell analysis and plotting of timetraces
clusterF = find(XI == 6);
ROI_cluster = numbers(:,idx(IX(clusterF)));
figure(93149);
for trx = 7:10
    subplot(1,4,trx-6); hold on;
    trial_timepoints = size(parameters.plane{1}.timetraces{trx},1);
    trial_timepoints2 = min(floor(trial_timepoints/7.5*100),size(parameters.pdgLUT{parameters.paradigms(trx)},1));
    cmap = distinguishable_colors(size(ROI_cluster,2));
    for k = 1:size(ROI_cluster,2)
        plot((1:trial_timepoints)/7.5+50/7.5-1.5,smooth(plane{ROI_cluster(1,k)}.timetraces{trx}(:,ROI_cluster(2,k))+k*100,5),'Color',cmap(k,:) );
    end
    xlabel('time [sec]'); ylabel('dF/F (arbitrary offset)')
    h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],100*[parameters.pdgLUT{parameters.paradigms(trx)}(1:10:trial_timepoints2,2); 0],'r');
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.4); hold on;
    h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],100*[parameters.pdgLUT{parameters.paradigms(trx)}(1:10:trial_timepoints2,3); 0],'b');
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.4);
    axis([min((1:trial_timepoints)/7.5+50/7.5-1.5) max((1:trial_timepoints)/7.5+50/7.5-1.5) -5 1395])
    grid on
end


%% show spatial distribution of a special cluster found above
clusterF = find(XI == 6);
ROI_cluster = numbers(:,idx(IX(clusterF)));

template = zeros(size(parameters.all_ROI));

parameters.planeColor = planeColor;

for k = 1:numel(clusterF)
    [x,y] = find(parameters.all_ROI == ROI_cluster(2,k) & parameters.planeColor ==ROI_cluster(1,k));
    for i = 1:numel(x)
        template(x(i),y(i)) = 1;
    end
end
figure(19); imagesc(template); axis equal off





%% select traces for hysteresis analysis and plot the traces
load('cluster_hysteresis_anaylsis.mat','clusterM')

for trx = 11:12
    traces = [];
    traces0 = [];
    numbers = [];
    qualitycheck = [];
    for k = 1:4
        x_length = size(plane{k}.timetraces_raw{trx},2);
        temp = plane{k}.timetraces{trx}(:,1:x_length);
        temp0 = plane{k}.timetraces{11}(:,1:x_length);
        traces = [traces, temp];
        traces0 = [traces0, temp0];
        retro_ix = [k*ones(x_length,1)'; 1:x_length];
        numbers = [numbers, retro_ix];
        qualitycheck = [qualitycheck, clusterM{k}];
    end
    %     idx = find(~isnan(sum(traces0)) & sum(traces0) ~= 0);
    idx = find(qualitycheck == 2);
    nb_clusters = 4;
    [~,XI,IX] = cluster_traces(traces0(:,idx),nb_clusters);
    traces_ordered = traces(:,idx(IX));

    for k = 1:size(traces_ordered,2)
        traces_ordered(:,k) = (traces_ordered(:,k) - mean(traces_ordered(:,k)))/std(traces_ordered(:,k));
    end
    % analysis of hysteresis effects
    figure(9999);
    subplot(1,2,trx-10);
    hold on;
    plot((1:size(traces_ordered,1))/7.5+50/7.5-1.5,conv2(traces_ordered,fspecial('gaussian',[3 1],1.5),'same') + repmat(linspace(0,90,size(traces_ordered,2)),[size(traces_ordered,1) 1]))
    xlabel('time [sec]'); ylabel('dF/F (arbitrary offset)')
    trial_timepoints = size(parameters.plane{1}.timetraces{trx},1);
    trial_timepoints2 = min(floor(trial_timepoints/7.5*100),size(parameters.pdgLUT{parameters.paradigms(trx)},1));
    h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[parameters.pdgLUT{parameters.paradigms(trx)}(1:10:trial_timepoints2,2); 0],'r');
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.4); hold on;
    h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[parameters.pdgLUT{parameters.paradigms(trx)}(1:10:trial_timepoints2,3); 0],'b');
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.4);
    axis([min((1:size(traces_ordered,1))/7.5+50/7.5-1.5) max((1:size(traces_ordered,1))/7.5+50/7.5-1.5) -5 95])
    grid on
end



%% pool and select interesting trials; plot correlations of state, maximum projection
% for each state (1D projection of the 2D corr matrix of states)

tracesZ = tracesX(1103:3606,:);
tracesZ((size(tracesZ,1)+1):(size(tracesZ,1)+851),:) = mean(cat(3,tracesX(3607:4457,:),tracesX(4458:5308,:)),3);
tracesZ((size(tracesZ,1)+1):(size(tracesZ,1)+851),:) = mean(cat(3,tracesX((3607:4457)+851*2,:),tracesX((4458:5308)+851*2,:)),3);
% tracesZ((size(tracesZ,1)+1):(size(tracesZ,1)+851),:) = tracesX((3607:4457)+851*2,:);
% tracesZ((size(tracesZ,1)+1):(size(tracesZ,1)+851),:) =tracesX((4458:5308)+851*2,:);

% normalization, if desired
tracesY = tracesZ;
% for j = 1:size(tracesY,2)
%     tracesY(:,j) = (tracesY(:,j) - mean(tracesY(:,j)))/std(tracesY(:,j));
% end
% use a short timetrace instead of instantaneous values; jj = 4 corresponds
% to a window of 9 frames, i.e., 1.2 sec
for jj = 4;% 10:50
    delayedTraces = tracesY';
    for k = 2:2:jj*2
        delayedTraces = [delayedTraces; circshift(tracesY',[k/2 k/2]); circshift(tracesY',[-k/2 -k/2]) ];
    end

    big_corr = corr(delayedTraces);

    % allin creates a sparse matrix that shows that maximum correlation
    % value in time
    small_corr = big_corr(1:2504,2505:end);
    allin = zeros(size(small_corr));
    for i = 1:size(small_corr,2)
        [max1(i),ix1] = max(small_corr(1:551,i));
        [max2(i),ix2] = max(small_corr(552:1102,i));
        [max3(i),ix3] = max(small_corr(1103:1803,i));
        [max4(i),ix4] = max(small_corr(1804:2504,i));
%         allin(ix,i) = 1;
    end
    figure(24); 
    subplot(2,2,1);
    plot((1:numel(max1))/7.5+50/7.5-1.5,max1,'r');hold on; % 1.5 sec is the delay between pump switch and odor arrival, empirically found for this experiment
    plot((1:numel(max1))/7.5+50/7.5-1.5,max2,'b');
    plot((1:numel(max1))/7.5+50/7.5-1.5,max3,'k');
    plot((1:numel(max1))/7.5+50/7.5-1.5,max4,'m');
    xlabel('time [sec]'); ylabel('Correlation value [0 1]')
    kkk = 7;
    axis([5.2 trial_timepoints/7.5 0 1]);  grid on
    subplot(2,2,3);
    trial_timepoints = size(parameters.plane{1}.timetraces{kkk},1);
    trial_timepoints2 = min(floor(trial_timepoints/7.5*100),size(parameters.pdgLUT{parameters.paradigms(kkk)},1));
    
    h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[0.3*parameters.pdgLUT{parameters.paradigms(kkk)}(1:10:trial_timepoints2,2); 0],'r');
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.4); hold on;
    h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[0.3*parameters.pdgLUT{parameters.paradigms(kkk)}(1:10:trial_timepoints2,3); 0],'b');
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.4);
    axis([5.2 trial_timepoints/7.5 0 1])
    hold off; grid on; xlabel('time [sec]'); ylabel('Correlation value [0 1]')
    subplot(2,2,2);
    plot((1:numel(max1))/7.5+50/7.5-1.5-851/7.5,max1,'r');hold on;
    plot((1:numel(max1))/7.5+50/7.5-1.5-851/7.5,max2,'b');
    plot((1:numel(max1))/7.5+50/7.5-1.5-851/7.5,max3,'k');
    plot((1:numel(max1))/7.5+50/7.5-1.5-851/7.5,max4,'m');
    axis([5.2 trial_timepoints/7.5 0 1]);  grid on
    xlabel('time [sec]'); ylabel('Correlation value [0 1]')
    kkk = 9;
    subplot(2,2,4);
    trial_timepoints = size(parameters.plane{1}.timetraces{kkk},1);
    trial_timepoints2 = min(floor(trial_timepoints/7.5*100),size(parameters.pdgLUT{parameters.paradigms(kkk)},1));
    
    h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[0.3*parameters.pdgLUT{parameters.paradigms(kkk)}(1:10:trial_timepoints2,2); 0],'r');
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.4); hold on;
    h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[0.3*parameters.pdgLUT{parameters.paradigms(kkk)}(1:10:trial_timepoints2,3); 0],'b');
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.4);
    axis([5.2 trial_timepoints/7.5 0.3 1])
    hold off; grid on; xlabel('time [sec]'); ylabel('Correlation value [0 1]')
end

% figure(9), imagesc(conv2(allin',[1 1 1; 1 1 1; 1 1 1],'same')); colormap(gray)
% figure(10), imagesc(small_corr'); colormap(gray)






%% pool and select interesting trials; plot correlations of state, maximum projection, Rainer's idea
% for each state (1D projection of the 2D corr matrix of states)

nb_frames = zeros(numel( plane{1}.timetraces),1);
for i = 1:numel(nb_frames)
    nb_frames(i) = size(plane{1}.timetraces{i},1);
end

test_trialx = [7 8];

for jjj = 1:numel(test_trialx)
    test_trial = test_trialx(jjj);
    
    template_trials = [3 4 5 6];
    template_windows = nb_frames(template_trials);
    trial_window = nb_frames(test_trial);

    tracesZ = [];
    for k = 1:numel(template_trials)
        tracesZ = [tracesZ; tracesX([(1:nb_frames(template_trials(k))) + sum(nb_frames(1:(template_trials(k)-1)))],:)];
    end
    X = zeros(size(tracesX([(1:nb_frames(test_trial(1))) + sum(nb_frames(1:(test_trial-1)))],:)));
    for k = 1:numel(test_trial)
        X = X + tracesX([(1:nb_frames(test_trial(k))) + sum(nb_frames(1:(test_trial(k)-1)))],:);
    end
    tracesZ = [tracesZ; X/numel(test_trial)];

    % normalization, if desired
    tracesY = tracesZ;
    if 0
        for j = 1:size(tracesY,2)
            tracesY(:,j) = (tracesY(:,j) - mean(tracesY(:,j)))/std(tracesY(:,j));
        end
    end
    idxs = find(~isnan(sum(tracesY)));
    tracesY = tracesY(:,idxs);

    % use a short timetrace instead of instantaneous values; jj = 4 corresponds
    % to a window of 9 frames, i.e., 1.2 sec
    delayedTraces = tracesY';
    %     for k = 2:2:jj*2
    %         delayedTraces = [delayedTraces; circshift(tracesY',[k/2 k/2]); circshift(tracesY',[-k/2 -k/2]) ];
    %     end

    % figure; hold on; plot(mean(delayedTraces),'g');  plot(std(delayedTraces)); 
    % [x, ~] = ginput(16);
    % x_selectors = round(x);
    % save('x_selectors.mat','x_selectors')

    big_corr = corr(delayedTraces);

    oxy{1} = x_selectors(1:4); % odor 1, block
    oxy{2} = x_selectors(5:8); % odor 2, block
    oxy{3} = x_selectors(9:12); % odor 1, gradient
    oxy{4} = x_selectors(13:16); % odor 2, gradient

    for j = 1:numel(oxy)
        start(j,:) = mean(big_corr((sum(template_windows)+1):end,oxy{j}(1):oxy{j}(2)),2);
        later(j,:) = mean(big_corr((sum(template_windows)+1):end,oxy{j}(2):oxy{j}(3)),2);
        off(j,:) = mean(big_corr((sum(template_windows)+1):end,oxy{j}(3):oxy{j}(4)),2);
    end

    % red=Phe, blue=Trp, black=Phe_v, magenta=Trp_v
    ccmap = [       1         0    0
        0         0         1
             0    0         0
             1         0    1];
    timet = (1:size(start,2))/7.5+50/7.5-1.5;
    figure(78+test_trialx(1)),
    subplot(1,3,1); hold on; for k = 1:4; plot(timet,smooth(start(k,:),5),'Color',ccmap(k,:)); end;
    xlabel('time [sec]'); ylabel('Correlation value [0 1]'); hold off; axis([min(timet) max(timet) -0.5 1])
    subplot(1,3,2); hold on; for k = 1:4; plot(timet,smooth(later(k,:),5),'Color',ccmap(k,:)); end;
    xlabel('time [sec]'); ylabel('Correlation value [0 1]'); hold off; axis([min(timet) max(timet) -0.5 1])
    subplot(1,3,3); hold on; for k = 1:4; plot(timet,smooth(off(k,:),5),'Color',ccmap(k,:)); end;
    xlabel('time [sec]'); ylabel('Correlation value [0 1]'); hold off; axis([min(timet) max(timet) -0.5 1])
end



%% show matrix of correlation of states
figure(7549), imagesc((0:4206)/7.5,(0:4206)/7.5,corr(tracesY'));
xlabel('time [sec]'); ylabel('time [sec]'); axis equal; axis([0 4206/7.5 0 4206/7.5]);


%% do a PCA and plot some of the components
[coeff,score] = pca(tracesZ);
figure(991); 
plot(score(:,1:6) + repmat(linspace(0,2000,6),[size(score,1) 1]))



