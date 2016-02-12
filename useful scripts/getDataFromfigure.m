% get data from figure object which is selected

h = gco;

a = get(h);
H{4,4} = a.CData;


% cmap = jet(999);
% cmap(500,:) = [ 1 1 1];
% figure(111);
% for pp = 1:4
%     D = cat(3,H{pp,1},H{pp,2},H{pp,3});
% %     figure(117),
% %     subplot(2,2,pp); hist(D(:),100);
%     R = max(max(D(:)),R);
%     P = min(min(D(:)),P);
%     
%     D = D + 0.6620;
%     D = D/(R-P+0.002)*2.5 - 0.77;
%     subplot(2,2,pp); imagesc( max(min(D,1),0)); axis off equal; %colormap(cmap)
% end

