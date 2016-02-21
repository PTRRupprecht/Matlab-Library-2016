

%% load mp3s

FolderList = dir();

for jj = 3:numel(FolderList)
    
    cd(FolderList(jj).name)

    FileList = dir('*.mp3');

    bitelength = 15750; % 2.5 sec ; one chunck will be 10 sec

    chunk_counter = 0;
    for i = 1:numel(FileList)
        disp(strcat('mp3 #',32,num2str(i),32,'out of',32,num2str(numel(FileList)),32,'within folder "',FolderList(jj).name,'".'));
        [y,Fs] = audioread(FileList(i).name);
        yy = zeros(floor(size(y,1)/7),2);
        for k = 1:7
            yy = yy + y((7:7:end)-k+1,:);
        end

        num_chunks = floor(size(yy,1)/bitelength)-3;

        chunked_song = zeros(bitelength*4,num_chunks*2);

        for p = 1:num_chunks
            chunked_song(:,(2*p-1):(2*p)) = yy((1:bitelength*4)+bitelength*(p-1),:);
        end
        
        chunk_counter = chunk_counter + num_chunks;
        
        save(strcat(FileList(i).name(1:end-4),'_chunks.mat'),'chunked_song');
    end
    
    cd ..
    Chunk_Counter(jj-2).filename = FolderList(jj).name;
    Chunk_Counter(jj-2).chun_counter = chunk_counter;
    
end

save('metaData.mat','Chunk_Counter')

%% read and play the chunks

Fileist = dir('*.mat');

for i = 1:50
    k = randi(numel(Fileist),1);
    F = load(Fileist(k).name);
    k = randi(size(F.chunked_song,2),1);
    soundsc(F.chunked_song(:,k),44100/7);
end


