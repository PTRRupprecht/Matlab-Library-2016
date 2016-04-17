
%% load data (plane, morphing paradigms, Filenames)
MatFileList = dir('Extracted*.mat');
load(MatFileList(1).name);
load('paradigms.mat');
F = dir('Fish1_*.tif');
FileNames = cat(1,F.name);


%% load paradigms
pdgLUT{1} = paradigm.p91; pdgLUT{2} = paradigm.p92;
pdgLUT{3} = paradigm.p93; pdgLUT{4} = paradigm.p94;
pdgLUT{5} = paradigm.p95; pdgLUT{6} = paradigm.p96;
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
figure(918), imagesc(RGBimage)
akZoom_PR();
set(gcf, 'WindowKeyPressFcn', {@OB_morphing_selectROI_9x,parameters});


