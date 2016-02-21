

clear all
% load list of files
% cd distorted
FileList = dir('AVG*.tif');

for k = 1:numel(FileList)
    clear meta
%     [A,result,meta.framerate,meta.zstep,meta.zoom,meta.motorpositions,meta.scalingfactors] = read_metadata_function(FileList(k).name);
    L{k} = imfinfo(FileList(k).name);
    meta.height = L{k}(1).Height;
    meta.width = L{k}(1).Width;
    meta.numberframes = numel(L{k});

    movie = read_movie(FileList(k).name,meta.width,meta.height,meta.numberframes,1,1,L{k},1);
    
    movie = bidi_align(movie);
%     movie = bidi_align_manual(movie,1);
%     figure(2), imagesc(movieX(:,:,15)); colormap(gray)
    movie = unwarp_precision(movie);
    
%     [px,py] = pixelsize_xy(meta.zoom,meta.width,meta.height);
    px = 0; py = px;
    
    filename = strcat(FileList(k).name(1:end-4),'_px',num2str(px*1000,'%4.0f'),'_py',num2str(py*1000,'%4.0f'),'.tif');
    imwrite(uint16(movie(:,:,1)),filename);
    for i = 2:size(movie,3)
        imwrite(uint16(movie(:,:,i)),filename,'WriteMode','append');
    end
end