

MatFileList = dir('Extracted*.mat');
load(MatFileList(1).name);

pdgLUT{1} = paradigm10Hz.p301; pdgLUT{2} = paradigm10Hz.p302;
pdgLUT{3} = paradigm10Hz.p32; pdgLUT{4} = paradigm10Hz.p32;
pdgLUT{5} = paradigm10Hz.p42; pdgLUT{6} = paradigm10Hz.p42;
nb_trials = size(plane{1}.anatomy,3);
smoothing = 5;
% find neuronal indizes
for j = 1:4; n_ix{j} = find(~isnan(sum(plane{j}.timetraces_raw{1}))); end

n_ix{2}

all_anatomy = zeros(1024);
all_anatomy(1:512,1:512) = plane{1}.anatomy(:,:,1);
all_anatomy(513:1024,1:512) = plane{2}.anatomy(:,:,1);
all_anatomy(1:512,513:1024) = plane{3}.anatomy(:,:,1);
all_anatomy(513:1024,513:1024) = plane{4}.anatomy(:,:,1);

all_anatomy(all_anatomy> quantile(all_anatomy(:),0.995)) = quantile(all_anatomy(:),0.995);




figure(918), imagesc(all_anatomy)



plane{1}.anatomy(:,:,1)



scaling = 2.5;
ROI_map_thresh = ROI_map; ROI_map_thresh(ROI_map_thresh>1) = 1;
RGBimage(:,:,1) = AVG_X*scaling; RGBimage(:,:,2) = AVG_X*scaling; RGBimage(:,:,3) = AVG_X*scaling;
RGBimage(:,:,2) = RGBimage(:,:,2) + ROI_map_thresh/scaling; RGBimage(:,:,1) = RGBimage(:,:,1) + ROI_map_thresh/scaling;
maxVal = max(max(RGBimage(:,:,1))); RGBimage(:,:,1) = RGBimage(:,:,1)/maxVal; RGBimage(:,:,2) = RGBimage(:,:,2)/max(max(RGBimage(:,:,2))); RGBimage(:,:,3) = RGBimage(:,:,3)/maxVal; 

set(handle1,'CData',RGBimage); caxis auto





plane_nb = 1;
neuron_nb = 2;

figure(41);
for k = 1:nb_trials
    % plot paradigm
    timestamps = (1:1900)/10;
    ax(k) = subplot(3,nb_trials,k);
    h = fill(timestamps,pdgLUT{k}(1:1900,1),'r');
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.4); hold on;
    h = fill(timestamps,pdgLUT{k}(1:1900,2),'b');
    set(h, 'EdgeColor','none', 'FaceAlpha', 0.4);
    axis([0 190 0 3])
    hold off;
    % plot response
    ax(k+nb_trials) = subplot(3,nb_trials,k+nb_trials); plot((1:1425)/7.5+50/7.5,smooth(plane{plane_nb}.timetraces{k}(1:1425,neuron_nb),smoothing),'k');
    axis([0 190 -10 205])
    
    % plot ROIs / neuron anatomy
    windowsize = 30;
    ROIX = squeeze(plane{plane_nb}.ROI_map(1,:,:));
    [x,y] = find(ROIX == neuron_nb); x = round(mean(x)); y = round(mean(y));
    xxx = max(1,x-windowsize):min(512,x+windowsize); yyy = max(1,y-windowsize):min(477,y+windowsize);
%     ROIXX = ROIX(xxx,yyy);
    ax(k+2*nb_trials) = subplot(3,nb_trials,k+2*nb_trials); imagesc(plane{plane_nb}.anatomy(xxx,yyy,k),[-30 70]); colormap(gray)
    title(FileNames(k).name,'Interpreter','None');
%     drawnow;
end




