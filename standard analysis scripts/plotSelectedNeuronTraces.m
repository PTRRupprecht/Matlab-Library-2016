function plotSelectedNeuronTraces(timetracesX_X,indizes,invers)

nbx = ceil(sqrt(length(indizes)));
nby = ceil(length(indizes)/nbx);

figure(41);
cmap = distinguishable_colors(length(indizes));
for i = 1:length(indizes)
    subplot(nbx,nby,i); hold on;
    if invers
        tat = numel(timetracesX_X):-1:1;
    else
        tat = 1:numel(timetracesX_X);
    end
    for k = tat
        plot(0.1:0.1:100,smooth(timetracesX_X{k}(:,indizes(i)),5),'Color',cmap(k,:));
    end
    box on;
    hold off;
end

end