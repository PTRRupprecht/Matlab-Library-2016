

jj = 1;
clear AVG
FileZ = dir('Fish1_OB2_pos3*.tif');
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




for pp = 1:4
    filename = strcat('AVGs_Fish1_OB2_pos3_',num2str(pp),'.tif');
    imwrite(uint16(squeeze(AVG(:,:,pp,1,jj))),filename);
    for kkk = 2:size(AVG,4)
        imwrite(uint16(squeeze(AVG(:,:,pp,kkk,jj))),filename,'WriteMode','append');
    end
end
