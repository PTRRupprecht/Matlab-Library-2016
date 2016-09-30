% Peter Rupprecht (c) 2016
% under MIT license

% list of folders
FolderList = dir();

counter_file = 1;
for jj = 3:9 % here, I had 7 folders with different composers
    
    cd(FolderList(jj).name)
    % list of mp3s
    FileList = dir('*.mp3');

    bitelength = 6300; % 1.0 sec; one chunck will be 10 sec

    chunk_counter = 0;
    for i = 1:numel(FileList)
        disp(strcat('mp3 #',32,num2str(i),32,'out of',32, ...
            num2str(numel(FileList)),32,'within folder "', ...
            FolderList(jj).name,'".'));
        % read mp3
        [y,Fs] = audioread(FileList(i).name);
        % downsample by factor 7
        yy = zeros(floor(size(y,1)/7),2);
        for k = 1:7
            yy = yy + y((7:7:end)-k+1,:);
        end
        yy = yy/7;

        % expected number of 10 sec junks
        num_chunks = floor(size(yy,1)/bitelength)-9;

        % generate junked pieces
        chunked_song = zeros(bitelength*10,num_chunks,'int16');
        for p = 1:num_chunks
            piece = yy((1:bitelength*10)+bitelength*(p-1),1);
            % 16bit, therefore multiply by 2^15-1=32767 for signed integer
            piece = piece*32767;
            chunked_song(:,p) = int16(piece);
        end
        % keep track of number of chunks for this folder
        chunk_counter = chunk_counter + num_chunks;
        % save as mat file (can be read in Python and Matlab)
        save(strcat(FileList(i).name(1:end-4),'_chunks.mat'),'chunked_song');
    end
    % go back 
    cd ..
    % write a mat file with metadata about the number of chunks contained
    % in the respective folder
    Chunk_Counter(counter_file).filename = FolderList(jj).name;
    Chunk_Counter(counter_file).chun_counter = chunk_counter;
    counter_file = counter_file + 1;
end

save('metaData.mat','Chunk_Counter')

%% read and play some of the the chunks
% move to folder with mat files first
FileList = dir('*.mat');
figure,
for i = 1:60
    k = randi(numel(FileList),1);
    F = load(FileList(k).name);
    k = randi(size(F.chunked_song,2),1);
    ding = double(F.chunked_song(:,k));
    ding = (ding - mean(ding))/std(ding);
    plot(ding+5*i);
    soundsc(ding,44100/7);
end


