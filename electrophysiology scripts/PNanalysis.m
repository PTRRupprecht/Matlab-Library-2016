

close all
%% show electrophysiological recordings
FileList = dir('*13*.xsg');
try; close 41; close 51; end
offsetV = 0;
offsetI = 0;
suppress_hum = 1;
cmap = lines(numel(FileList));
counter = 1;
VCindex = [];
CCindex = [];
for i = 1
    load(FileList(i).name,'-mat');
    A = data.ephys.trace_1;
    A2 = A;
    if 0 % subtract 50 Hz noise
        window = 10000; % 2 sec
        for kk = 1:numel(A)/window;
            X = A((1:window) + (kk-1)*window);
            X = X - mean(X);
            A((1:window) + (kk-1)*window) = A((1:window) + (kk-1)*window) - repmat(mean( reshape(X,[200 numel(X)/200]),2),[numel(X)/200  1]);
        end
    end
    samplerate = header.ephys.ephys.sampleRate;
    timet = (1:numel(A))/samplerate;%*1000; % ms
   
    
    transients = zeros(1500,400);
    for k = 1:400
        transients(:,k) = A((1:1500)+(k-1)*1500);
    end
    
end

template10 = mean(transients(:,1:100),2) - mean(mean(transients(1:240,1:100),2));
template70 = mean(transients(:,101:200),2) - mean(mean(transients(1:240,101:200),2));
template10x = mean(transients(:,201:300),2) - mean(mean(transients(1:240,201:300),2));
template35 = mean(transients(:,301:400),2) - mean(mean(transients(1:240,301:400),2));

% raw data
figure(51), plot(template10); hold on; plot(template10x);
plot(template70,'r'); plot(template35,'k')
% corrected data
figure(52), plot(template70-template10*7,'r'); hold on; plot(template35 - template10x*3.5,'k');

% XXX
timecourse2D = transients(:,101:200) - 0*repmat(mean(transients(:,1:100),2)*7,[1 100]);
figure(63);  hold on; cmap = jet(100);
for k = 1:100
    plot(smooth(timecourse2D(:,k),5),'Color',cmap(k,:));
end
figure(64); plot(mean(transients(1000:1200,:),1))


% corrected data, timecourse - +70 mV
timecourse2D = transients(:,101:200) - repmat(mean(transients(:,1:100),2)*7,[1 100]);
figure(53); subplot(1,2,1); hold on; cmap = jet(100);
for k = 1:100
    plot(smooth(timecourse2D(:,k),1),'Color',cmap(k,:));
end
subplot(1,2,2); imagesc(timecourse2D(250:1250,:))
% corrected data, timecourse - +35 mV
timecourse2D = transients(:,301:400) - repmat(mean(transients(:,201:300),2)*3.5,[1 100]);
figure(54); subplot(1,2,1); hold on; cmap = jet(100);
for k = 1:100
    plot(smooth(timecourse2D(:,k),5),'Color',cmap(k,:));
end
subplot(1,2,2); imagesc(timecourse2D(250:1250,:))
% corrected data, timecourse - 0 mV
timecourse2D = transients(:,201:300) - repmat(mean(transients(:,001:100),2),[1 100]);
figure(55); subplot(1,2,1); hold on; cmap = jet(100);
for k = 1:100
    plot(smooth(timecourse2D(:,k),5),'Color',cmap(k,:));
end
subplot(1,2,2); imagesc(timecourse2D(250:1250,:))


