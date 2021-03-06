clear all
global clut2b timetracesX_X ROI_map_X movie_AVG_X
load clut2b

%% load list of files
FileList_0 = dir('Fish1_odor2Lys_00*_.tif');
FileList = FileList_0;
clear meta
[A,result,meta.framerate,meta.zstep,meta.zoom,meta.motorpositions,meta.scalingfactors] = read_metadata_function(FileList_0(1).name);
for kkk = 1:numel(FileList)
    L{kkk} = imfinfo(FileList(kkk).name);
    meta.height = L{kkk}(1).Height;
    meta.width = L{kkk}(1).Width;
    meta.numberframes(kkk) = numel(L{kkk});
end


%% read raw data from hard disk
binning = 1;
meta.framerate = meta.framerate/binning;
nb_frames_total = sum(floor(meta.numberframes/binning));
clear movie movie_p counter_planes
useful_range_start = 200;


nb_planes = 4;
meta.framerate = meta.framerate/nb_planes;


nb_frames_perplane = floor(nb_frames_total/nb_planes);
for pp = 1:4 % which plane do you want to choose right now?
    counter_planes{pp} = 0;
    clear movie_p
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

    %% subtract preamp ringing
    template_window = [];
    for k = 1:8; template_window = [template_window; [1:25]'+400*(k-1)]; end
    template = mean(movie_p{pp}(:,:,template_window),3);
    template_odd = mean(template(1:2:end,:),1);
    template_even = mean(template(2:2:end,:),1);
    for k = 1:size(template,2)/2
        template(2*k-1,:) = template_odd;
        template(2*k,:) = template_even;
    end
    for k = 1:size(movie_p{pp},3)
        movie_p{pp}(:,:,k) = movie_p{pp}(:,:,k) - template;
    end

    % unwarp full stack
    movie_p{pp} = unwarp_precision(movie_p{pp});

    %% AVG images
    clear AVG
    for k = 1:8
        AVG(:,:,k) = mean(movie_p{pp}(:,:,(51:400)+(k-1)*400),3);
        result_conv =fftshift(real(ifft2(conj(fft2(AVG(:,:,k))).*fft2(AVG(:,:,1)))));
        [y,x] = find(result_conv==max(result_conv(:))); %Find the 255 peak
        result_conv =fftshift(real(ifft2(conj(fft2(AVG(:,:,1))).*fft2(AVG(:,:,1)))));
        [y0,x0] = find(result_conv==max(result_conv(:))); %Find the 255 peak
        offsety(k) = y-y0;
        offsetx(k) = x-x0;
    end


    %% treat each trial separately
    % odor switches at frame 100
    figure(854); hold on;
    for trial_nb = 1:8
        movie_trial = circshift(movie_p{pp}(:,:,(50:400)+(trial_nb-1)*400),[offsety(trial_nb) offsetx(trial_nb) 0]);
        AVG_movie(:,:,trial_nb) = mean(movie_trial,3);

        % calculate activity maps
        offset = -25;
        f0_window = [1 100];
        response_window = [118 150];
        plot1 = 0; plot2 = 0; DF_movie_yesno = 0; % figure number
        [DF_reponse(:,:,trial_nb),DF_master(:,:,trial_nb),DF_movie] = dFoverF(movie_trial,offset,meta.framerate,plot1,plot2,DF_movie_yesno,f0_window,response_window);
        % local correlation map (computational slightly expensive, but helpful)
        tilesize = 16;
        localCorrelations(:,:,trial_nb) =zeros(size(DF_master(:,:,trial_nb)));
    end

    %% final ROIs, selected by comparing with all 8 trials in imageJ
    clear timetracesXX timetracesXX_raw
    for trial_nb = 1:8
        offset = -25;
        movie_trial = circshift(movie_p{pp}(:,:,(50:400)+(trial_nb-1)*400),[offsety(trial_nb) offsetx(trial_nb) 0]);
        df_scale = [-20 100];
        AVG_Z = AVG_movie(:,:,trial_nb);
        AVG_Z(AVG_Z>80) = 80;
        ROI_map_input = squeeze(plane{pp}.ROI_map(trial_nb,:,:));
        [ROI_mapX(trial_nb,:,:),timetracesX,timetracesX_raw] = timetraces_singleplane_recovery(movie_trial,AVG_Z,offset,DF_reponse(:,:,trial_nb),DF_master(:,:,trial_nb),localCorrelations(:,:,trial_nb),df_scale,ROI_map_input,meta,1,AVG_Z);
        timetracesXX(:,:,trial_nb) = timetracesX;
        timetracesXX_raw(:,:,trial_nb) = timetracesX_raw;
    end

    plane{pp}.ROI_map = ROI_mapX;
    plane{pp}.timetraces = timetracesXX;
    plane{pp}.timetraces_raw = timetracesXX_raw;
    plane{pp}.DF_reponse = DF_reponse;
    plane{pp}.meta = meta;
    plane{pp}.anatomy = AVG_movie;

end
save(strcat('Extracted_Data_fish_01-03-16_Lys.mat'),'plane');

