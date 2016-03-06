
% Plotting overview of neuronal activity
% - one subplot per neuron with 8 repetitions
% - one figure per experiment

MatFileList = dir('Fish*X.mat');
% for all 7 experiments
for yy = 1:numel(MatFileList)
    load(MatFileList(yy).name);
    clear fluoXT fluoX
    % for all 8 trials
    for trailx = 1:size(plane{1}.timetraces,3)
        fluoX = [];
        planesIX = [];  
        planesIXX = [];
        % for all 4 planes
        for jj = 1:numel(plane)
            fluotrace = plane{jj}.timetraces_raw;
            indizes = find(~isnan(sum(fluotrace(:,:,1),1)));
            fluotrace = fluotrace(:,indizes,:);
            for kk = 1:size(fluotrace,2)
                ftrace = smooth(fluotrace(:,kk,trailx),25); F0 = min(ftrace(25:end-25));
                fluotraceX = (fluotrace(:,kk,trailx) - F0)/F0*100;
                fluoX = [fluoX, smooth(fluotraceX,5)];
                planesIX = [planesIX,jj];
                planesIXX = [planesIXX,indizes(kk)];
            end
        end
        fluoXT(:,:,trailx) = fluoX;
    end
    % plot each experiment separately; each neuron gets a subplot
    figure(1218+yy); hold on; title(MatFileList(yy).name)
    for trailx = 1:size(fluoXT,2)
        subplot(ceil(size(fluoXT,2)/floor(sqrt(size(fluoXT,2)))),floor(sqrt(size(fluoXT,2))),trailx); imagesc((1:351)/7.5,[],squeeze(fluoXT(:,trailx,:))',[0 200]);
        title(strcat('Plane',32,num2str(planesIX(trailx)),', neuron',32,num2str(planesIXX(trailx))));
    end
    xlabel(strcat('time (sec)',32,32,32,32,32,32,32,MatFileList(yy).name),'Interpreter','None'); ylabel('trial nb');
end

