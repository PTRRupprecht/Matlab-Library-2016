
FileList = dir('*.xsg');

for i = 1% :numel(FileList)

    load(FileList(i).name,'-mat');
    A = data.ephys.trace_1;
    samplerate = header.ephys.ephys.sampleRate;
    timet = (1:numel(A))/samplerate*1000; % ms
end
figure(1), plot(timet, smooth(A,10))