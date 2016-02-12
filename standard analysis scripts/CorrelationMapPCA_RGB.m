
%% PCA and locate PCs spatially

[COEFF,SCORE] = princomp(fullX);

cmap = jet(999);
figure(811);
D = repmat(ROI_mapXX,[1 1 3]);
for j = 1:max(D(:))
    [x,y] = find(ROI_mapXX == j);
    if ~isempty(x)
        A1 = corr(timetracesX_X_raw(:,j),SCORE(:,1));
        A2 = corr(timetracesX_X_raw(:,j),SCORE(:,2));
        A3 = corr(timetracesX_X_raw(:,j),SCORE(:,3));
        for k = 1:numel(x)
            D(x(k),y(k),:) = cat(3,A1,A2,A3);
        end
    end
end
D = (D+1)/2; D(D == 0.5) = [0];
imagesc(abs(D),[0 1]); axis off equal; %colormap(cmap)

