
cd('M:\rupppete\electrophysiology2016\');

FileList = dir('*.xsg');

for i = 1:numel(FileList)
    i
    load(FileList(i).name,'-mat');
    A = data.ephys.trace_1;
    samplerate = header.ephys.ephys.sampleRate;
    timet = (1:numel(A))/samplerate;%*1000; % ms
    trace = smooth(A,10);
    figure(1), plot(timet(1:10:end), trace(1:10:end))
    pause(0.5);
end
step = samplerate*150;
A(end+1:step*5) = 0;
for j = 1:5
    junkedTrace(j,:) = A((1:step)+(j-1)*step);
end
figure(413);hold on; cmap = lines(10);
for j = 1:5
    plot((1:step)/samplerate,smooth(junkedTrace(j,:)+20*j,50),'Color',cmap(j,:));
%     pause(0.5);
end
plot((1:150*1e2)/1e2,waveform4(1:150*1e2,2)*38,'k');


figure(1), plot(timet, smooth(A,10))