

clear all
% load list of files
% cd distorted
FileList = dir('*Fish1_GAD_12dpf_TectumCoarse_zstack_001_*.tif');

for k = 1:numel(FileList)
    clear meta
%     [A,result,meta.framerate,meta.zstep,meta.zoom,meta.motorpositions,meta.scalingfactors] = read_metadata_function(FileList(k).name);
    L{k} = imfinfo(FileList(k).name);
    meta.height = L{k}(1).Height;
    meta.width = L{k}(1).Width;
    meta.numberframes = numel(L{k});

    movie = read_movie(FileList(k).name,meta.width,meta.height,meta.numberframes,1,1,L{k},1);
    [~,~,~,zstep,zoom,~,~] = read_metadata_function(FileList(k).name);
    
    %% subtract preamp ringing
    if 1
        template_window = [];
        for kk = 1:numel(FileList); template_window = 1:meta.numberframes; end
        template = mean(movie(:,:,template_window),3);
        template_odd = mean(template(1:2:end,:),1);
        template_even = mean(template(2:2:end,:),1);
        for kk = 1:size(template,1)/2
            template(2*kk-1,:) = template_odd;
            template(2*kk,:) = template_even;
        end
        template = template - imfilter(template,fspecial('gaussian',[15 15],7),'replicate');
        % figure, imagesc(template)
        for kk = 1:size(movie,3)
            movie(:,:,kk) = movie(:,:,kk) - template;
        end
    end
    
    movie = bidi_align(movie);
%     movie = bidi_align_manual(movie,1);
%     figure(2), imagesc(movie(:,:,15)); colormap(gray)
    movie = unwarp_precision(movie);
    
%     figure(1); imagesc(movie(:,:,80));caxis([8400 8490]); colormap(gray)
%     figure(2); imagesc(unwarp_precision(movie(:,:,80))); caxis([8400 8490]); colormap(gray)
    
    [px,py] = pixelsize_xy(zoom,meta.width,meta.height);
%     px = 0; py = px;
    
    filename = strcat(FileList(k).name(1:end-4),'_pixelsizey_',num2str(py*1000,'%4.0f'),'_zstepsize_',num2str(zstep*1000,'%6.0f'),'.tif');
    imwrite(uint16(movie(:,:,1)),filename);
    for i = 2:size(movie,3)
        imwrite(uint16(movie(:,:,i)),filename,'WriteMode','append');
    end
end
toc