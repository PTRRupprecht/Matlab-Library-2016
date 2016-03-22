
% select neuron to be viewed across trials for closer inspection

function OB_morphing_selectROI(~,event,parameters)
    clear keyword
    keyword = event.Key;
    switch(keyword)
        case 'x'
            [x,y] = ginput(1); x = round(x); y = round(y);
            plane_nb = parameters.planeColor(y,x);
            neuron_nb = parameters.all_ROI(y,x);
            figure(41);
            for k = 1:parameters.nb_trials
                % plot paradigm
                ax(k) = subplot(3,parameters.nb_trials,k);
                h = fill(parameters.timestamps,parameters.pdgLUT{k}(1:1900,1),'r');
                set(h, 'EdgeColor','none', 'FaceAlpha', 0.4); hold on;
                h = fill(parameters.timestamps,parameters.pdgLUT{k}(1:1900,2),'b');
                set(h, 'EdgeColor','none', 'FaceAlpha', 0.4);
                axis([0 190 0 3])
                plot((1:1425)/7.5+50/7.5,smooth(parameters.plane{plane_nb}.timetraces{k}(1:1425,neuron_nb),parameters.smoothing)/70,'k');

                hold off;
                % plot response
                ax(k+parameters.nb_trials) = subplot(3,parameters.nb_trials,k+parameters.nb_trials); plot((1:1425)/7.5+50/7.5,smooth(parameters.plane{plane_nb}.timetraces{k}(1:1425,neuron_nb),parameters.smoothing),'k');
                axis([0 190 -10 205])

                % plot ROIs / neuron anatomy
                windowsize = 30;
                ROIX = squeeze(parameters.plane{plane_nb}.ROI_map(1,:,:));
                [x,y] = find(ROIX == neuron_nb); x = round(mean(x)); y = round(mean(y));
                xxx = max(1,x-windowsize):min(512,x+windowsize); yyy = max(1,y-windowsize):min(477,y+windowsize);
                ax(k+2*parameters.nb_trials) = subplot(3,parameters.nb_trials,k+2*parameters.nb_trials); imagesc(parameters.plane{plane_nb}.anatomy(xxx,yyy,k),[-30 70]); colormap(gray)
                title(parameters.FileNames(k).name,'Interpreter','None');
            end
            figure(918),
    end
end