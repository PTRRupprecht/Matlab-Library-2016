% INPUT is a cell with each of the n cells consisting of T x N matrizes, with T the
% timepoints for each trial n and N the number of neurons
function [corrMap,timepoints] = corrNaNs(INPUT)

    clear timepoints;
    for j = 1:numel(INPUT)
        timepoints(j) = size(INPUT{j},1);
    end
    corrMap = zeros(size(cell2mat(INPUT'),1));
    for k = 1:numel(INPUT)
        for j = 1:numel(INPUT)
            indizesx = find(~isnan(sum(INPUT{k},1)));
            indizesy = find(~isnan(sum(INPUT{j},1)));
            indizes = intersect(indizesy,indizesx);
            corrMap( (sum(timepoints(1:k-1))+1):sum(timepoints(1:k)),(sum(timepoints(1:j-1))+1):sum(timepoints(1:j)) ) = corr(INPUT{k}(:,indizes)',INPUT{j}(:,indizes)');
        end
    end
    
end