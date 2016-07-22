


%% show electrophysiological recordings
FileList = dir('*14*.xsg');
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
   
    
    AA = zeros(1500,4);
    for k = 1:4
        for j = 1:100
            AA(:,k) = A((1:1500)+(j-1)*1500 + 150000*(k-1));
        end
        AA(:,k) = AA(:,k)/j;
        AA(:,k) = AA(:,k) - mean(AA(1:240,k));
    end
    figure, plot(smooth(-AA(:,3)+AA(:,4)/3.5,7))
    hold on; plot(smooth(-AA(:,1)+AA(:,2)/7,7),'k')
    
    FF(:,i-21) = AA - mean(AA);
    
    % end
    %     
    % WW = mean(FF(:,1:2),2); WW2 = mean(FF(:,3:4),2)/5; WW3 = mean(FF(:,5:6),2)/8;
    % WW = smooth(WW,5);
    % WW2 = smooth(WW2,5);
    % WW3 = smooth(WW3,5);
    % figure, hold on; plot(WW); plot(WW2+2,'r'); plot(WW3+4,'k')
    % hold on; plot(WW2-WW+6,'r'); plot(WW3-WW+8,'k')

 
 
    trace = smooth(AA,5);
    figure(431), subplot(4,4,i-20); plot(AA)
    timet = (1:numel(AA))/samplerate;%*1000; % ms
    
    counter = counter + 1;
    if ~strcmp(header.ephys.ephys.amplifierSettings.Amp_700B_1.mode,'I-Clamp')
        figure(41), hold on; plot(timet(1:10:end), trace(1:10:end)+offsetV,'Color',cmap(i,:))
        text(timet(end)+1,median(trace(1:10:end)+offsetV),num2str(i),'FontSize',12)
        offsetV = offsetV + max(trace)-min(trace);
        VCindex = [VCindex; i];
    else
        figure(51), hold on; plot(timet(1:10:end), trace(1:10:end)+offsetI,'Color',cmap(i,:))
        text(timet(end),median(trace(1:10:end)+offsetI),num2str(i),'FontSize',12)
        offsetI = offsetI + max(trace)-min(trace);
        CCindex = [CCindex; i];
    end
end
VCindex = VCindex'
CCindex = CCindex'
