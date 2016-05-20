function [traces_ordered,XI,IX] = cluster_traces(traces,nb_clusters)

    FF = traces;
    dF = FF;
    corrMatrix = corr(dF);
    corrMatrix = corrMatrix - eye(size(corrMatrix));
    dissimilarity = 1 - corrMatrix';

    % perform complete linkage clustering
    Z = linkage(dissimilarity,'complete');
    if nb_clusters == 0
        IX = optimalleaforder(Z,round((dissimilarity-eye(size(dissimilarity)))*1000)/1000);
        XI = ones(size(IX));
    else
        % decide on a cutoff
        cutoff = 10;
        l_groups = 1;
        while l_groups < nb_clusters
            cutoff = cutoff - 0.1;
            groups = cluster(Z,'cutoff',cutoff,'criterion','distance');
            l_groups = length(unique(groups));
        end
        [XI,IX] = sort(groups);
    end
    traces_ordered = dF(:,IX);
end