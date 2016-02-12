%% usage
% global clusterX
% for k = 1:size(timetracesX_X_raw,2)
%     handle1 = figure(2198); plot(smooth(timetracesX_X_raw(:,k),5));
%     set(gcf, 'units','normalized','Position',[0.1 0.1 0.55 0.620]);
%     set(gcf, 'WindowKeyPressFcn', {@chooseCluster,k,handle1});
%     grid on;
%     waitfor(gcf);
% end
% save('ClusterX.mat','ppClusterX')


function chooseCluster(~,event,kk,handle1)

    global clusterX;
    clear keyword
    keyword = event.Key;
    switch(keyword)
        case 'leftarrow'
            clusterX(kk) = 1;
        case 'rightarrow'
            clusterX(kk) = 2;
        case 'uparrow'
            clusterX(kk) = 3;
        case 'downarrow'
            clusterX(kk) = 4;
    end
    disp(num2str(kk));
    close 2198;
end