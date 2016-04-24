function [traces_ordered,XI,IX] = cluster_traces(traces,nb_clusters)

    FF = traces;
    dF = FF;
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
    traces_ordered = dF(:,IX);
end