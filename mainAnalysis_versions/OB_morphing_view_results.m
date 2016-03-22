
%% load data (plane, morphing paradigms, Filenames)
MatFileList = dir('Fish1_extr*.mat');
load(MatFileList(1).name);

%% load paradigms
pdgLUT{1} = paradigm10Hz.p301; pdgLUT{2} = paradigm10Hz.p302;
pdgLUT{3} = paradigm10Hz.p32; pdgLUT{4} = paradigm10Hz.p32;
pdgLUT{5} = paradigm10Hz.p42; pdgLUT{6} = paradigm10Hz.p42;
nb_trials = size(plane{1}.anatomy,3);
smoothing = 5;
% find neuronal indizes
% for j = 1:4; n_ix{j} = find(~isnan(sum(plane{j}.timetraces_raw{1}))); end

%% make anatomy clickable (RBG overlay image)
all_anatomy = zeros(1024);
all_anatomy(1:512,1:512) = plane{1}.anatomy(:,:,1); all_ROI(1:512,1:512) = squeeze(plane{1}.ROI_map(1,:,:));
all_anatomy(513:1024,1:512) = plane{2}.anatomy(:,:,1); all_ROI(513:1024,1:512) = squeeze(plane{2}.ROI_map(1,:,:));
all_anatomy(1:512,513:1024) = plane{3}.anatomy(:,:,1); all_ROI(1:512,513:1024) = squeeze(plane{3}.ROI_map(1,:,:));
all_anatomy(513:1024,513:1024) = plane{4}.anatomy(:,:,1); all_ROI(513:1024,513:1024) = squeeze(plane{4}.ROI_map(1,:,:));
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
parameters.timestamps = (1:1900)/10;
parameters.pdgLUT = pdgLUT;
parameters.plane = plane;
parameters.FileNames = FileNames;
parameters.all_ROI = all_ROI;
parameters.planeColor = planeColor;
parameters.smoothing = smoothing;

%% figure and callback function
figure(918), imagesc(RGBimage)
akZoom_PR();
set(gcf, 'WindowKeyPressFcn', {@OB_morphing_selectROI,parameters});


