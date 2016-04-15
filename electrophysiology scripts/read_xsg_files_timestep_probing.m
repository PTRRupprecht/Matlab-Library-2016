
cd('M:\rupppete\electrophysiology2016\');

FileList = dir('*.xsg');

for i = 16
    i
    load(FileList(i).name,'-mat');
    A = data.ephys.trace_1;
    samplerate = header.ephys.ephys.sampleRate;
    timet = (1:numel(A))/samplerate;%*1000; % ms
    trace = smooth(A,30);
    figure(1), plot(timet(1:10:end), trace(1:10:end))
    pause(1.5);
end
step = samplerate;
% A(end+1:step*5) = 0;
for j = 1:60
    junkedTrace(j,:) = A((1:step)+(j-1)*step+7500);
end
figure(8), imagesc([1:step]/10,[],junkedTrace); colorbar; xlabel('time [ms]');

y = mean(junkedTrace,1);
x = [1:step]/10;
cftool % use Levenherdt-Marquardt
figure(9); plot([1:step]/10,mean(junkedTrace,1))



