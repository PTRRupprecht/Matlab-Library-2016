clear all
global clut2b timetracesX_X ROI_map_X movie_AVG_X
load clut2b

%% load list of files
FileList_0 = dir('Fish1_10_pdgm96_001*.tif');
FileList = FileList_0;
clear meta
[A,result,meta.framerate,meta.zstep,meta.zoom,meta.motorpositions,meta.scalingfactors] = read_metadata_function(FileList_0(1).name);
for kkk = 1:numel(FileList)
    L{kkk} = imfinfo(FileList(kkk).name);
    meta.height = L{kkk}(1).Height;
    meta.width = L{kkk}(1).Width;
    meta.numberframes(kkk) = numel(L{kkk});
end

% check whether acquisitions have been aborted or not
nb_planes = 4;
for i = 1:numel(meta.numberframes)
    if rem(meta.numberframes(i),nb_planes) ~= 0
         meta.numberframes(i) = floor(meta.numberframes(i)/nb_planes)*nb_planes;
    end
end


%% read raw data from hard disk
binning = 1;
meta.framerate = meta.framerate/binning;
meta.framerate = meta.framerate/nb_planes;
nb_frames_total = sum(floor(meta.numberframes/binning));

clear movie movie_p counter_planes
nb_frames_perplane = floor(nb_frames_total/nb_planes);
for pp = 1:4 % which plane do you want to choose right now?
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
end


movie_pp = zeros(2*size(movie_p{1},1),2*size(movie_p{1},2),size(movie_p{1},3));

movie_pp(1:512,1:512,:) = movie_p{1};
movie_pp(1:512,(1:512)+512,:) = movie_p{3};
movie_pp((1:512)+512,1:512,:) = movie_p{2};
movie_pp((1:512)+512,(1:512)+512,:) = movie_p{4};


%% subtract preamp ringing
template_window = [];
for k = 1:numel(FileList); template_window = [template_window; [1:25]'+sum(meta.numberframes(1:(k-1)))/nb_planes]; end
template = mean(movie_pp(:,:,template_window),3);
template_odd = mean(template(1:2:end,:),1);
template_even = mean(template(2:2:end,:),1);
for k = 1:size(template,2)/2
    template(2*k-1,:) = template_odd;
    template(2*k,:) = template_even;
end
% figure, imagesc(template)
for k = 1:size(movie_pp,3)
    movie_pp(:,:,k) = movie_pp(:,:,k) - template;
end

% unwarp full stack
% movie_p{pp} = unwarp_precision(movie_p{pp});

%% AVG images
clear AVG
for trial_nb = 1:numel(FileList)
    AVG(:,:,trial_nb) = mean(movie_pp(:,:,(51:meta.numberframes(trial_nb)/nb_planes)+sum(meta.numberframes(1:(trial_nb-1)))/nb_planes),3);
    result_conv =fftshift(real(ifft2(conj(fft2(AVG(:,:,trial_nb))).*fft2(AVG(:,:,1)))));
    [y,x] = find(result_conv==max(result_conv(:))); %Find the 255 peak
    result_conv =fftshift(real(ifft2(conj(fft2(AVG(:,:,1))).*fft2(AVG(:,:,1)))));
    [y0,x0] = find(result_conv==max(result_conv(:))); %Find the 255 peak
    offsety(trial_nb) = y-y0;
    offsetx(trial_nb) = x-x0;
end


%% treat each trial separately
% odor switches at frame 100
figure(854); hold on;
for trial_nb = 1:numel(FileList)
    movie_trial = movie_pp(:,:,(50:meta.numberframes(trial_nb)/nb_planes)+sum(meta.numberframes(1:(trial_nb-1)))/nb_planes);
    AVG_movie(:,:,trial_nb) = mean(movie_trial,3);

    % calculate activity maps
    offset = -30;
    f0_window = [1 100];
    response_window = [300 min(550,size(movie_trial,3))];
    plot1 = 0; plot2 = 0; DF_movie_yesno = 0; % figure number
    [DF_reponse(:,:,trial_nb),DF_master(:,:,trial_nb),DF_movie] = dFoverF(movie_trial,offset,meta.framerate,plot1,plot2,DF_movie_yesno,f0_window,response_window);
    % local correlation map (computational slightly expensive, but helpful)
    tilesize = 16;
%     localCorrelations(:,:,trial_nb) = zeros(size(DF_reponse(:,:,trial_nb)));
    localCorrelations(:,:,trial_nb) = localCorrelationMap(movie_trial,tilesize);

    subplot(2,ceil(numel(FileList)/2),trial_nb); imagesc(DF_reponse(:,:,trial_nb),[-0.5 2])
%     subplot(2,ceil(numel(FileList)/2),trial_nb); imagesc(AVG_movie(:,:,trial_nb),[-30 70]); colormap(gray)
end
akZoom('all_linked')

%% preliminary ROIs got from test trials
ROI_map_input = zeros(size(AVG_movie(:,:,1)));
trial_nb = 1; 
offset = -25;
movie_trial = circshift(movie_pp(:,:,(50:meta.numberframes(trial_nb)/nb_planes)+sum(meta.numberframes(1:(trial_nb-1)))/nb_planes),[offsety(trial_nb) offsetx(trial_nb) 0]);
df_scale = [-20 100];
AVG_Z = AVG_movie(:,:,trial_nb);
AVG_Z(AVG_Z>80) = 80;

[ROI_mapX(trial_nb,:,:),timetracesX,timetracesX_raw] = timetraces_singleplane(movie_trial,AVG_Z,offset,DF_reponse(:,:,trial_nb),DF_master(:,:,trial_nb),localCorrelations(:,:,trial_nb),df_scale,ROI_map_input,meta,1,AVG_Z);
% figure, imagesc(conv2(timetracesX,fspecial('gaussian',[25 1],23),'same'));
ROI_map_input = squeeze(ROI_mapX(trial_nb,:,:));


%% calculate correlation map based on a time window and a time course from one of the ROIs
window = 300:700;
AA = timetracesX_raw(window,1);
xcorr_tot = zeros(size(movie_trial,1),size(movie_trial,2));
for j = 1:size(movie_trial,1)
    if mod(j,50) == 0; disp(j/size(movie_trial,1)); end
    xcorr_tot(j,:) = corr(AA,squeeze(movie_trial(j,:,window))');
end

figure, imagesc(xcorr_tot); colormap(gray)






