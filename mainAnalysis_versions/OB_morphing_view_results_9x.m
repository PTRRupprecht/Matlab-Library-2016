
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
select_trial = 5;
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


% trial
tracesX = [];
for trx = 1:13
    x_length = size(plane{1}.timetraces{7},2);
    traces = [];
    for k = 1:4
        temp = plane{k}.timetraces{trx}(:,1:x_length);
        temp0 = plane{k}.timetraces{8}(:,1:x_length);
        idx = find(~isnan(sum(temp0)) & sum(temp0) ~= 0);
        traces = [traces, temp(:,idx)];
    end
figure(2)
    subplot(4,4,trx), imagesc((1:size(traces,1))/7.5+50/7.5,[],traces',[-40 400]);
%     subplot(4,6,trx*2),imagesc(corr(traces'));
    
figure(4)
    trial_timepoints = size(parameters.plane{1}.timetraces{trx},1);
    trial_timepoints2 = min(floor(trial_timepoints/7.5*100),size(parameters.pdgLUT{parameters.paradigms(trx)},1));
    subplot(4,4,trx),  h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[parameters.pdgLUT{parameters.paradigms(trx)}(1:10:trial_timepoints2,2); 0],'r');
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.4); hold on;
    h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[parameters.pdgLUT{parameters.paradigms(trx)}(1:10:trial_timepoints2,3); 0],'b');
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.4);

   
end


