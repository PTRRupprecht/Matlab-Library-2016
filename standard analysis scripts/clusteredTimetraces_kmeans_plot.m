%% detect clusters (k-means) and plot the ordered timetraces


select_traces = find(clusterX == 2);

FF = fullX;

nb_clusters = 7; % choose number of clusters to detect

detect_traces = find(nansum(FF));
dF = FF(:,detect_traces);
corrMatrix = corr(dF);
corrMatrixOld = corrMatrix;
corrMatrix = corrMatrix - eye(size(corrMatrix));
dissimilarity = 1 - corrMatrix';

% perform complete linkage clustering
Z = linkage(dissimilarity,'complete');
% decide on a cutoff
cutoff = 10;
l_groups = 1;
while l_groups < nb_clusters
    cutoff = cutoff - 0.1;
    groups = cluster(Z,'cutoff',cutoff,'criterion','distance');
    l_groups = length(unique(groups));
end

[XI,IX] = sort(groups);
dF_ordered = dF(:,IX);

corrMatrix2 = corr(dF_ordered);

% figure, imagesc(dF_ordered)
% [IX,XI] = sort(max(dF_ordered));

cmap = distinguishable_colors(size(dF_ordered,2));
figure(4); hold on;
for j = 1:size(dF_ordered,2)
    plot(smooth(dF_ordered(:,IX(end-j-1)) + 440*j,4),'Color',cmap(j,:));
end
    