

FileX = dir('1603*_gilad');

for jj = 1:numel(FileX)
    cd(FileX(jj).name);
    FileZ = dir('*pdg*.tif');
    for kkk = 1:numel(FileZ)
        jj
        kkk
        L{kkk} = imfinfo(FileZ(kkk).name);
        clear X
    	X = read_movieLX(FileZ(kkk).name,512,512,numel(L{kkk}),1,1,L{kkk},1);
        for pp = 1:4
            AVG(:,:,pp,kkk,jj) = mean(X(:,:,(200:4:end-5)+pp),3);
        end
    end
    cd ..
end


for jj = 1 3%:numel(FileX)
    for pp = 1:4
        filename = strcat('FishI_Phe_Trp_plane',num2str(pp),'.tif');
        imwrite(uint16(squeeze(AVG(:,:,pp,1,jj))),filename);
        for kkk = 2:6
            imwrite(uint16(squeeze(AVG(:,:,pp,kkk,jj))),filename,'WriteMode','append');
        end
    end
end
