

%% load list of files
filename = 'Fish2_FoodPos2_001_.tif';
FileList_0 = dir(filename);
FileList = FileList_0;
clear meta
[AAA,result,meta.framerate,meta.zstep,meta.zoom,meta.motorpositions,meta.scalingfactors] = read_metadata_function(FileList_0(1).name);
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

for pp = 1:4 
%% subtract preamp ringing
template_window = [];
for k = 1:numel(FileList); template_window = [template_window; [1:size(movie_p{pp},3)]'+sum(meta.numberframes(1:(k-1)))/nb_planes]; end
template = mean(movie_p{pp}(:,:,template_window),3);
template_odd = mean(template(1:2:end,:),1) - mean(template(1:1:end,:),1);
template_even = mean(template(2:2:end,:),1) - mean(template(1:1:end,:),1);
for k = 1:size(template,2)/2
    template(2*k-1,:) = template_odd;
    template(2*k,:) = template_even;
end
% figure, imagesc(template)
for k = 1:size(movie_p{pp},3)
    movie_p{pp}(:,:,k) = movie_p{pp}(:,:,k) - template;
end

% unwarp full stack
movie_p{pp} = unwarp_precision(movie_p{pp});

reference = mean(movie_p{pp}(:,:,200:300),3);

[movie_p{pp},offsety_resolved,offsetx_resolved] = alginWithinTrial(movie_p{pp},reference);

end

for pp = 1:4
    bytesPerPixel = 8; % for double

    headerString = AAA(1).ImageDescription;
    F = cellstr(headerString);
    ts = TifStream(strcat(filename(1:end-4),num2str(pp),'_forRF','.tif'),size(movie_p{pp},2),size(movie_p{pp},1),F{1});

    % read & write; ca. 8-10 ms for a 512x512 frame
    for k = 25:500%size(movie,3)
        A=movie_p{pp}(:,:,k);
        ts.appendFrame(A);

    end
    % close pipes
    ts.close();
end


