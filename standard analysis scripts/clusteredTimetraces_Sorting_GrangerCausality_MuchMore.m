%% read in data

Matrix = Sheet;
start_index = 10;
offset = 10;

%% calculate dF/F

dF = zeros(6000,63);
for i = 1:63
    timetrace = Matrix(:,(i-1)*4+3) - 10;
    center(1,i) = Matrix(1,(i-1)*4+4);
    center(2,i) = Matrix(1,(i-1)*4+5);
    F0 = median(timetrace(start_index:end));
    dF(:,i) = smooth((timetrace-F0)/F0,10)*100;
end
figure(652), imagesc(dF'); colorbar

%% cluster analysis
background = zeros(512);
nb_clusters = 3;

[idx,Centroid,~,~] = kmeans(dF',nb_clusters);
cmap = lines(nb_clusters);
figure(413), imagesc(background); colormap(gray); hold on;
for k = 1:63
    plot(center(1,k),center(2,k),'.','Color',cmap(idx(k),:),'MarkerSize',38)
end
figure(6849), hold on;
for k = 1:nb_clusters
    plot(Centroid(k,:)+30*(k-1),'Color',cmap(k,:))
end


%% optimal leaf ordering for clusters
cmap = distinguishable_colors(5);
counter = 0;
colorcounter = [1 3 4 5 2];
smoothing = 3;
figure(95); hold on;
for j = [1:5]
    jj = j;
    good_indizes = find(xy_index == j & clusterX == 2);
    if j == 5;
        good_indizes = [good_indizes, find(xy_index == j+1 & clusterX == 2)];
    end
    traces{jj} = timetracesX_cooked(:,good_indizes);
    
    corrMatrix = corr(traces{jj});
    corrMatrixOld = corrMatrix;
    corrMatrix = corrMatrix - eye(size(corrMatrix));
    dissimilarity = 1 - corrMatrix';

    Z = linkage(dissimilarity,'complete');

    leafOrder = optimalleaforder(Z,round((dissimilarity-eye(size(dissimilarity)))*1000)/1000);

    for i = 1:numel(good_indizes)
        plot((1:size(timetracesX_cooked,1))/meta.framerate,smooth(170*counter+timetracesX_cooked(:,good_indizes(leafOrder(i))),smoothing),'Color',cmap(colorcounter(j),:));
        counter = counter + 1;
        fulltraces(:,counter) = timetracesX_cooked(:,good_indizes(leafOrder(i)));
    end
end


%% correlation matrix and cluster tree

corrMatrix = corr(dF);

figure, imagesc(corrMatrix); hold on;
for k = 1:63
    plot(k,k,'.','Color',cmap(idx(k),:),'MarkerSize',36)
end

%% additional cluster analysis & sorting of the corr-matrix based on these correlations

%# remove diagonal elements
corrMatrix = corrMatrix - eye(size(corrMatrix));
%# and convert to a vector (as pdist)
dissimilarity = 1 - corrMatrix';

%# decide on a cutoff
%# remember that 0.4 corresponds to corr of 0.6!
cutoff = 0.5; 

%# perform complete linkage clustering
Z = linkage(dissimilarity,'complete');

cutoff = 4;
figure(4111); dendrogram(Z,0,'colorthreshold',cutoff)

%# group the data into clusters
cmap = lines(63);
nb_groups = 0;
breakout_nb = 2;
for i = 1:1000
    cutoff = (1-i/1000)*max(Z(:,3));
    groups = cluster(Z,'cutoff',cutoff,'criterion','distance');
    if nb_groups ~= length(unique(groups));
        figure(413), imagesc(background); colormap(gray); hold on;
        for k = 1:63
            plot(center(1,k),center(2,k),'.','Color',cmap(groups(k),:),'MarkerSize',38)
        end
    end
    nb_groups = length(unique(groups));
    if breakout_nb == length(unique(groups))
        break;
    end
end

[A,IX] = sort(groups);
dF_ordered = dF(:,IX);

corrMatrix2 = corr(dF_ordered);
figure, imagesc(corrMatrix2); hold on;
for k = 1:63
    plot(k,k,'.','Color',cmap(A(k),:),'MarkerSize',36)
end

figure(413), imagesc(background); colormap(gray); hold on;
for k = 1:63
    plot(center(1,IX(k)),center(2,IX(k)),'.','Color',cmap(A(k),:),'MarkerSize',38)
end


%% distance in space correlates with functional correlation ?
distMap = squareform(pdist(center(:,IX)'));
figure, imagesc(distMap)
figure, scatter(distMap(:),1-corrMatrix2(:))

%% time-dependent clustering (state reservoir)

corrMatrix_states = corr(dF');
figure(54), imagesc(corrMatrix_states);

%% time-dependent correlation matrix
figure,
for i = 1:6
    corrMatrix = corr(dF_ordered((i-1)*1000+1:i*1000,:));
    subplot(2,3,i); imagesc(corrMatrix,[-0.2 0.5]); drawnow; pause(1);
end


%% firing rate estimation

firing = dF;
level_cut = 10;

firing(firing<level_cut) = 0;
firing(firing>0) = 1;

firing_idx = sum(firing);

cmmap = gray(max(firing_idx));
figure(413), imagesc(background); colormap(gray); hold on;
for k = 1:63
    plot(center(1,IX(k)),center(2,IX(k)),'.','Color',cmmap(firing_idx(k),:),'MarkerSize',38)
end


%% mutual information

nb_bins = 25; % nb of activity levels
AA = dF;
binlocations = zeros(nb_bins,size(AA,2)); % left edge
P_marginals = zeros(nb_bins,size(AA,2)); % simple distributions
for i = 1:size(AA,2)
    [P_marginals(:,i), binlocations(:,i)] = hist(AA(:,i),nb_bins); 
end

P_marginals = P_marginals + 0.001;

PXY = zeros(nb_bins,nb_bins,size(AA,2),size(AA,2));
for i = 1:size(AA,2)
    i
    for j = 1:size(AA,2) % 1:size(AA,2) % 
        for k = 1:nb_bins
            for m = 1:nb_bins
                if k < nb_bins
                    [posk, ~] = find(AA(:,i) >= binlocations(k,i) & AA(:,i) < binlocations(k+1,i));
                else
                    [posk, ~] = find(AA(:,i) >= binlocations(k,i));
                end
                if m < nb_bins
                    [posm, ~] = find(AA(:,j) >= binlocations(m,j) & AA(:,j) < binlocations(m+1,j));
                else
                    [posm, ~] = find(AA(:,j) >= binlocations(m,j));
                end
                PXY(i,j,k,m) = length(intersect(posk,posm));
            end
        end
    end
end

PXY = PXY + 0.001;

Mutual_Information = zeros(size(AA,2),size(AA,2));
for i = 1:size(AA,2)
    i
    for j = 1:size(AA,2)
        
        for k = 1:nb_bins
            for m = 1:nb_bins
                Mutual_Information(i,j) = Mutual_Information(i,j) + PXY(i,j,k,m)/1000*log2( PXY(i,j,k,m)/1000/ ( P_marginals(k,i)*P_marginals(m,j)/1000^2));
            end
        end
    end
end

MAP = tril(Mutual_Information,-2);

figure; imagesc(Mutual_Information,[min(MAP(:)), max(MAP(:))])


%% Granger causality

nlags = 4;
[ret] = cca_granger_regress(dF',nlags);

[PR] = cca_findsignificance(ret,0.01,1);
GC = GC.*PR;
cca_plotcausality(GC);







