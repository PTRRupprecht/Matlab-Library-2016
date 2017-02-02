
%% load dataset (dataset1), which had been created using the script create_dataset1.m
% The main purpose of this script is to bring together several experiments
% with sometimes different odor onset and framerate.
% The dataset consists of single trial elements, each associated with a
% experimental series ('dataset') and an odorant ('odor')

load('OB_tuning_enhanced.mat');

% find out experiments (= unique elements in the dataset)
for k = 1:numel(dataset1)
    setList{k} = dataset1{k}.dataset;
end
[datasetList, A,B] = unique(setList); 

c2ounter = 1;
c1ounter = 1;
for kkk = [1:5 7 9:numel(datasetList)] % discard experiments that do not show odor discrimination in the neuronal responses

    setIX = find(strcmp(datasetList{kkk},setList));
    % discard those experiments that do not have clear onset/offset
    % (e.g. morphing experiments)
    for j = 1:numel(setIX)
        if isnan(dataset1{setIX(j)}.onset)
            setIX(j) = 0;
        elseif isnan(dataset1{setIX(j)}.offset)
            setIX(j) = 0;
        end
    end
    setIX(setIX == 0) = [];
    
    % find maximum number of timepoints
    maxTime = min( size(dataset1{setIX(1)}.timetraces_raw{1},1) , size(dataset1{setIX(2)}.timetraces_raw{1},1) ); % shorten every timetrace to maximum
    
    clear XX_timetraces XX_odor
    for j = 1:numel(setIX) % go through the trials of an experiment
        X = dataset1{setIX(j)}; % load dataset
        X_timetraces = cell2mat(X.timetraces_raw); % concatenate timetraces
        if any(strcmp('cluster',fieldnames(X))) % any preselection of stable neurons; applies to morphing experiments mostly
            cluster = cell2mat(X.cluster);
            indizes = find(cluster ~= 2);
            X_timetraces(:,indizes) = [];
        end
        for kk = 1:size(X_timetraces,2) % calculate dF/F from the raw data
            ftrace = smooth(X_timetraces(:,kk),25); F0 = min(ftrace(25:end-25));
            X_timetraces(:,kk) = (X_timetraces(:,kk) - F0)/F0*100;
        end
        XX_timetraces(:,:,j) = X_timetraces(1:maxTime,:); % save everything in a bigger matrix
        XX_odor{j} = X.odor; % save odor identity
    end
    
    % calculate correlations between trials across the diagonal, similar to
    % distance of trial trajectories in PCA space
    clear X1 X2 CorrVec
    counter = 1;
    for j = 1:(numel(setIX)-1)
        for k = (j+1):numel(setIX)
            indxs = find(~isnan(sum(XX_timetraces(:,:,j))) & ~isnan(sum(XX_timetraces(:,:,k))));
            CorrMat = corr(XX_timetraces(:,indxs,j)',XX_timetraces(:,indxs,k)');
            CorrMat = conv2(CorrMat,fspecial('gaussian',[5 5],dataset1{setIX(1)}.meta.framerate),'same');
            CorrVec(:,counter) = diag(CorrMat);
            X1{counter} = XX_odor{j}{1};
            X2{counter} = XX_odor{k}{1};
            counter = counter + 1;
        end
    end
    nb_neurons(kkk) = size(XX_timetraces(:,indxs,j),2);
    % save everything in a super-structure
    CorrVector{kkk}.X1 = X1;
    CorrVector{kkk}.X2 = X2;
    CorrVector{kkk}.CorrVec = CorrVec;
    CorrVector{kkk}.onset = dataset1{setIX(1)}.onset;
    CorrVector{kkk}.offset = dataset1{setIX(1)}.offset;
    CorrVector{kkk}.framerate = dataset1{setIX(1)}.meta.framerate;
    CorrVector{kkk}.timetraces = XX_timetraces;
    
    % onset determintation based on refinement of the assumed onset
    on = round(CorrVector{kkk}.onset*CorrVector{kkk}.framerate) - 25;
    dd = nanmean(nanmean(CorrVector{kkk}.timetraces(on:(on+50),:,:),2),3);
    thresh = (max(dd)-min(dd))/2+min(dd);
    correction = find(dd > thresh,1,'first');
    CorrVector{kkk}.onset = CorrVector{kkk}.onset + (correction - 25)/CorrVector{kkk}.framerate;
    CorrVector{kkk}.offset = CorrVector{kkk}.offset + (correction - 25)/CorrVector{kkk}.framerate;
    
    % plot single traces of each trial pair, with color encoding if the two
    % trials have used the same odorant or not
    timet = (1:maxTime)/CorrVector{kkk}.framerate - CorrVector{kkk}.onset;
    figure(31); subplot(4,4,kkk); hold on; 
    for mm = 1:size(CorrVec,2)
        if strcmp(X1(mm),X2(mm))
            plot(timet,CorrVec(:,mm),'k') % same odorant
        else
            plot(timet,CorrVec(:,mm),'r') % different odorant
        end
    end
    axis([-10 50 -0.2 1])
end


%% plot correlation matrix of all and everything


AllTraces = [];
for kkk = [5 7 9:numel(datasetList)]
    timet = (1:size(CorrVector{kkk}.timetraces,1))/CorrVector{kkk}.framerate - CorrVector{kkk}.onset;
    
    shiftX = round((CorrVector{kkk}.onset-10)*CorrVector{kkk}.framerate);
    ttimetraces = circshift(CorrVector{kkk}.timetraces,[-shiftX 0 0]);
    
    
    if kkk < 5
        ttimetraces(:,:,5:6) = NaN;
    end
    clear ttimetracesJ
    if kkk == 5
        for j = 1:6
            ttimetracesJ(:,:,j) = resample(ttimetraces(:,:,j),2,1);
        end
    else
        ttimetracesJ = ttimetraces;
    end


    AllTraces = [AllTraces,ttimetracesJ(1:500,:,:)];

end

% PCA plot
AllTraces = reshape(permute(AllTraces,[1 3 2]),[500*6 467]);
ixs = ~isnan(sum(AllTraces(:,:,:)));
AllTracesY = AllTraces(:,ixs);

[coeff,score,latent] = pca(AllTracesY);

figure, plot(score(:,1:6))

colors = {'r', 'k', 'b', 'r', 'k', 'b'};
figure, hold on;
kk = [1 2 3];
for i = 1:6
    timepoints = (50:330)+(i-1)*500;
    plot3(smooth(score(timepoints,kk(1)),3),smooth(score(timepoints,kk(2)),3),smooth(score(timepoints,kk(3)),3),'Color',colors{i});
end


% use AllTraces without the modifications of the last paragraph
AllTracesX = reshape(AllTraces,[500 467*6]);
ixs = find(~isnan(sum(sum(AllTracesX(:,:,:)),3)));
AllTracesX = AllTracesX(:,ixs,:);


CorrMatt = corr(squeeze(AllTracesX(:,:))');

[x,y] = find(CorrMatt> 0.58 & CorrMatt < 0.62);


txime = ((1:numel(CorrMatt(j,1:end))))/7.5 - 10;
figure, imagesc(txime,txime,corr(squeeze(AllTracesX(:,:))'),[-0.2 1])
axis equal; axis([min(txime) max(txime) min(txime) max(txime)]);
colormap(brewermap(256,'Greys'))
cmap = jet(55);
hold on; 
figure, contourf(txime,txime,CorrMatt,[ -.2 0.8])



%% plot nice anatomy map

figure(4132),
for j = 1:12
    for h = 1:4
        subplot(4,12,j + 12*(h-1)); imagesc(dataset1{A(j)}.anatomy{h},[-3 45]); % axis equal off;
        colormap(gray)
    end
end







%% average over all correlation vectors (decorrelation / distance measurement)

sameX = [];
diffX = [];
weightD = [];
weightS = [];
counterS = 1;
counterD = 1;
for k = [1:5 7 9:numel(datasetList)] %1:numel(CorrVector) % go through experiments
    clear W % W will be the matrix containing the correlation trace of a trial pair
    if k == 5 % special treatment for recordings with lower framerate
        H = CorrVector{k}.CorrVec;
        W(2:2:2*size(H,1),:) = H;
        W((2:2:2*size(H,1))+1,:) = H;
    else
        W = CorrVector{k}.CorrVec;
    end
    weight(k) = size(CorrVector{k}.timetraces,2); % weight = number of neurons in this trial
    onset = round(-CorrVector{k}.onset*7.5+27*7.5);
    sameXX = NaN*ones(951,120); % final structure
    diffXX = NaN*ones(951,120);
    counterSS = 1; counterDD = 1;
    for jj = 1:numel(CorrVector{k}.X1) % add correlation vectors with appropriate temporal shift
        if strcmp(CorrVector{k}.X1(jj),CorrVector{k}.X2(jj))
            sameXX(onset:(onset+size(W,1)-1),counterSS) = W(:,jj);
            counterSS = counterSS + 1;
        else
            diffXX(onset:(onset+size(W,1)-1),counterDD) = W(:,jj);
            counterDD = counterDD + 1;
        end
    end
    mean(nanmean(sameXX(200:400,2)) - nanmean(diffXX(200:400,2))) % if this value is very low, the recorded neurons cannot properly distinguish odors
    if nanmean(nanmean(sameXX(200:400,2)) - nanmean(diffXX(200:400,2))) > 0.20 % arbitrary threshold, does not discard anything here
        diffX = [diffX, diffXX];
        weightD = [weightD, ones(size(diffXX,2),1)*weight(k)];
        sameX = [sameX, sameXX];
        weightS = [weightS, ones(size(sameXX,2),1)*weight(k)];
    end
end
% eliminate NaNs
id = find(~isnan(diffX(300,:))); diffX = diffX(:,id); weightD = weightD(id);
is = find(~isnan(sameX(300,:))); sameX = sameX(:,is); weightS = weightS(is);

% weighted mean of all correlation traces
timeF = (1:951)/7.5-27;
DmeanTrace = sum(diffX.*repmat(weightD,[size(diffX,1) 1]),2)/sum(weightD);
SmeanTrace = sum(sameX.*repmat(weightS,[size(sameX,1) 1]),2)/sum(weightS);
DstdTrace = sqrt(sum((diffX-repmat(DmeanTrace,[1 size(diffX,2)])).^2.*repmat(weightD,[size(diffX,1) 1]),2)/sum(weightD));
SstdTrace = sqrt(sum((sameX-repmat(SmeanTrace,[1 size(sameX,2)])).^2.*repmat(weightS,[size(sameX,1) 1]),2)/sum(weightS));

% plot correlation trace for same/different odors
figure(25), subplot(1,2,1); plot(timeF,sum(diffX.*repmat(weightD,[size(diffX,1) 1]),2)/sum(weightD),'k')
axis([-10 50 0.2 0.5])
subplot(1,2,2);  plot(timeF,sum(sameX.*repmat(weightS,[size(sameX,1) 1]),2)/sum(weightS),'k');
axis([-10 50 0.2 0.8])

% figure(54); shadedErrorBar(timeF,DmeanTrace,DstdTrace,'k',1) % uses
% toolbox/fileExchange (Rob Cambell)

figure(99), imagesc(diffX)

%% plot timetraces

TT = CorrVector{9}.timetraces;
TT = cat(2,TT,CorrVector{10}.timetraces);
TT = cat(2,TT,CorrVector{11}.timetraces);
TT = cat(2,TT,CorrVector{12}.timetraces);
ix = find(~isnan(sum(sum(TT(:,:,:),3),1)));
[~, IX, XI] = cluster_traces(squeeze(TT(:,ix,1)),7);
TTT = TT(:,ix(XI),:);
figure(7), imagesc([TTT(:,:,1)' TTT(:,:,2)' TTT(:,:,3)' TTT(:,:,4)' TTT(:,:,5)' TTT(:,:,6)' ]);

TI = [];
for j = 9:12
TT = CorrVector{j}.timetraces;
ix = find(~isnan(sum(sum(TT(:,:,:),3),1)));
[~, IX, XI] = cluster_traces(squeeze(TT(:,ix,1)),7);
TTT = TT(:,ix(XI),:);
TI = [TI, TTT];

figure(j), imagesc([TTT(:,:,1)' TTT(:,:,2)' TTT(:,:,3)' TTT(:,:,4)' TTT(:,:,5)' TTT(:,:,6)' ]);
end

figure(41),
for k = 1:6
    [~,~,ixx] = cluster_traces([TI(:,:,1); TI(:,:,2); TI(:,:,3); TI(:,:,4); TI(:,:,5); TI(:,:,6)],6);
    subplot(1,6,k); plot(timet,conv2(TI(:,ixx,k),fspecial('gaussian',[6 1],4),'same') + (ones(1,size(TI,1))'*linspace(1,10000,size(TI,2))) )
    axis([-10 58 -50 10200])
end

%% plot correlation of states thing
clear StateXX
counter = 1;
weightsum = 0;
for kkk = [1:5 7 9:numel(datasetList)]
    if kkk == 5 % upsample 8 planes recording
        clear W
        H = CorrVector{kkk}.timetraces;
        W(2:2:2*size(H,1),:,:) = H;
        W((2:2:2*size(H,1))+1,:,:) = H;
        TX = W;
    else
        TX = CorrVector{kkk}.timetraces;
    end
    if 0 % normalize timetraces
        for j = 1:size(TX,2)
            for i = 1:size(TX,3)
                TX(:,j,i) = (TX(:,j,i) - mean(TX(:,j,i)))/std(TX(:,j,i));
            end
        end
    end
    
    for jj = 1:size(TX,3) % calculate correlation matrix and align temporally
        G = TX(:,:,jj);
        ix = find(~isnan(sum(G)));
        StateX = corr(G(:,ix)');
        StateX = circshift(StateX,10*7.5-round([CorrVector{kkk}.onset*7.5 CorrVector{kkk}.onset*7.5]));
        
        StateXX(:,:,counter) = weight(kkk)*StateX(1:503,1:503);
        weightsum = weightsum + weight(kkk);
        counter = counter + 1;
    end
end
avgStates = sum(StateXX,3)/weightsum;
% clear autoCorr
% autoCorr = zeros(300);
% for j = 1:300
%     autoCorr(j,end-j+1:end) = avgStates(j,1:j);
% end
% autoCorr(autoCorr==0) = NaN;
% figure, plot(autoCorr(1:20:end,:)')

figure(42), % plot averaged (auto-)correlation matrix
imagesc([1:300]/7.5-10,[1:300]/7.5-10,avgStates(1:300,1:300),[-0.2 1])


%% single neuron tuning over time

clear DiffIndex DiffTempX
for kkk = [1:5 7 9:numel(datasetList)]
    if kkk == 5 % upsample 8 planes recording
        clear W
        H = CorrVector{kkk}.timetraces;
        W(2:2:2*size(H,1),:,:) = H;
        W((2:2:2*size(H,1))+1,:,:) = H;
        TX = W;
    else
        TX = CorrVector{kkk}.timetraces;
    end
    
    % re-extract odor sequence ...
    waiter = 0; odorsequence = {};
    counterU = 1;
    counterStart = numel(CorrVector{kkk}.X1) + 1;
    odorsequence{1} = CorrVector{kkk}.X2{end};
    while waiter == 0
        try
            counterStart = counterStart - counterU;
            odorsequence{end+1} = CorrVector{kkk}.X1{counterStart};
            counterU = counterU + 1;
        catch
            waiter = 1;
        end
    end
    odorsequence = fliplr(odorsequence);
    clear DiffTemp
    % calculate for each neuron (jj) and each timepoint (tt) the average
    % (over mm and nn, the different trials) difference index between
    % responses to same and different stimuli (DiffIndex)
    for jj = 1:size(TX,2)
        TR = squeeze(TX(:,jj,:));
        for tt = 1:size(TR,1)
            TV = squeeze(TR(tt,:));
            countD = 0; Df = 0;
            countS = 0; Sa = 0;
            for mm = 1:(numel(TV)-1)
                for nn = (mm+1):numel(TV)
                    if strcmp(odorsequence{mm},odorsequence{nn})
                        countS = countS + 1;
                        Sa = Sa + abs(TV(mm) - TV(nn));
                    else
                        countD = countD + 1;
                        Df = Df + abs(TV(mm) - TV(nn));
                    end
                end
            end
            DiffTemp(jj,tt) = Df/countD - Sa/countS;
        end
    end
    
    % sort timetraces according to difference index during the first 5
    % seconds ('DiffTempX', to be used later for pooled data)
    onset = round(CorrVector{kkk}.onset*7.5);
    onset2 = onset + round(7.5*5);
    onset3 = 300;
    onset4 = 350;
    
    DiffTempX{kkk} = mean(DiffTemp(:,onset:onset2),2);
    if kkk == 5 
        DiffTempY{kkk} = mean(DiffTemp(:,(onset3:onset4)-50),2);
    else
        DiffTempY{kkk} = mean(DiffTemp(:,onset3:onset4),2);
    end    
    DiffIndex{kkk} = DiffTemp;
end

% plot ordered timetraces
timeTV = (1:503)/7.5 - onset/7.6;
% DDD = [DiffIndex{5}(:,1:503); DiffIndex{9}; DiffIndex{10}; DiffIndex{11}; DiffIndex{12}];
DDD = [circshift(DiffIndex{5}(:,1:503),[0 50]); DiffIndex{7}(:,1:503); DiffIndex{9}(:,1:503); DiffIndex{10}(:,1:503); DiffIndex{11}(:,1:503); DiffIndex{12}(:,1:503)];
% DDT = [DiffTempX{9}; DiffTempX{10}; DiffTempX{11}; DiffTempX{12}];
DDT = [DiffTempY{5}; DiffTempY{7}; DiffTempY{9}; DiffTempY{10}; DiffTempY{11}; DiffTempY{12}];
[~, IX] = sort(DDT);
figure(574), imagesc(timeTV,[],DDD(IX(1:353),:),[-40 400])
