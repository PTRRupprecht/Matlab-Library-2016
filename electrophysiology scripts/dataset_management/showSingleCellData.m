%% show annotation viewer
cd('C:\Data\rupppete\PhD\electrophysiology2016 - Copy\AnnotationViewer');
load('StacksUint16.mat');
figure(99);
imshow3DofDp(Horizontal,Sagittal,'NeuronList_Annotation',[8387 8557]);


%% show stacks for this cell
stacks = dir('*.tif');
try; close 1; close 2; close 3; close 4; close 5; close 6; end; drawnow;
for j = 1:numel(stacks)
    stack_info = imfinfo(stacks(j).name);
    [~,~,~,zstep,zoom,~,~] = read_metadata_function(stacks(j).name);
    [movie,movie_AVG] = read_movie(stacks(j).name,stack_info(1).Width,stack_info(1).Height,numel(stack_info),1,1,stack_info,1);
    figure(j); imshow3Dfull(movie,[8463-40 8463+40]);
    disp(['Stack nb.',32,num2str(j),32,'with z-step size',32,num2str(zstep),32,'um.']);
end

% show full Dp stack
backdoor = pwd; date = backdoor(end-10:end-5);
cd('C:\Data\rupppete\PhD\electrophysiology2016 - Copy\Dp_stacks');
listofstacks = dir(strcat(date,'*.tif'));
try; close 1111; end
try
    stack_info = imfinfo(listofstacks(1).name);
    [movie,movie_AVG] = read_movie(listofstacks(1).name,stack_info(1).Width,stack_info(1).Height,numel(stack_info),1,1,stack_info,1);
    [~,~,~,zstep,zoom,~,~] = read_metadata_function(listofstacks(1).name);
    figure(1111); imshow3Dfull(movie,[8463-40 8463+40])
    disp(['Overview stack with z-step size',32,num2str(zstep),32,'um.']);
catch
    disp('No full stack available');
end
cd(backdoor)


%% show electrophysiological recordings
FileList = dir('*.xsg');
try; close 41; close 51; end
offsetV = 0;
offsetI = 0;
suppress_hum = 1;
cmap = lines(numel(FileList));
counter = 1;
VCindex = [];
CCindex = [];
for i = 1:numel(FileList)
    load(FileList(i).name,'-mat');
    A = data.ephys.trace_1;
    A2 = A;
    if 1 % subtract 50 Hz noise
        window = 10000; % 2 sec
        for kk = 1:numel(A)/window;
            X = A((1:window) + (kk-1)*window);
            X = X - mean(X);
            A((1:window) + (kk-1)*window) = A((1:window) + (kk-1)*window) - repmat(mean( reshape(X,[200 numel(X)/200]),2),[numel(X)/200  1]);
        end
    end
    samplerate = header.ephys.ephys.sampleRate;
    timet = (1:numel(A))/samplerate;%*1000; % ms
    trace = smooth(A,10);
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

