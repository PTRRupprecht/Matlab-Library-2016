
%% PCA and locate PCs spatially, using ROIs and extracted timetraces

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

%% correlation maps using PCAs etc.

[COEFF,SCORE] = princomp(timetracesX(:,1:60));

x = 1:900;
AA(:,1) = zeros(900,1); AA(140:190,1) = 1;
AA(:,2) = zeros(900,1); AA(1:120,2) = 1;
AA(:,3) = zeros(900,1); AA(150:618,3) = 1; 
figure(99), plot(AA)
 
for k = 1:3
    k
    xcorr_tot = zeros(size(movie,1),size(movie,2));
    for j = 1:size(movie,1)
        if mod(j,50) == 0; disp(j/size(movie,1)); end
        xcorr_tot(j,:) = corr(AA(:,k),squeeze(movie(j,:,:))');
    end
    H{k} = xcorr_tot;
end
A = H{1}; B = H{2}; C = H{3};
D = cat(3,A,B,C);
D = (D+1)/2;
D = D*2.5; 
figure(31); imagesc(min(max(D-0.75,0),1)); axis off equal; %colormap(cmap)



figure(41); hold on;
for k = 1:9
    subplot(3,3,k); imagesc(H{k},[-0.3 0.3]); colormap(gray(512));
end

