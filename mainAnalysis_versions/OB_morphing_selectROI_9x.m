
% select neuron to be viewed across trials for closer inspection

function OB_morphing_selectROI_9x(~,event,parameters)
    clear keyword
    keyword = event.Key;
    switch(keyword)
        case 'x'
            [x,y] = ginput(1); x = round(x); y = round(y);
            plane_nb = parameters.planeColor(y,x);
            neuron_nb = parameters.all_ROI(y,x);
            figure(41);
            nb_y = 4;
            nb_x = 6;
            for k = 1:parameters.nb_trials
                % plot paradigm
                row_pos = floor(k/nb_x-0.001);
                trial_timepoints = size(parameters.plane{plane_nb}.timetraces{k},1);
                trial_timepoints2 = min(floor(trial_timepoints/7.5*100),size(parameters.pdgLUT{parameters.paradigms(k)},1));
                ax(k) = subplot(nb_y,nb_x,k + row_pos*nb_x);
                h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[parameters.pdgLUT{parameters.paradigms(k)}(1:10:trial_timepoints2,2); 0],'r');
                set(h, 'EdgeColor','none', 'FaceAlpha', 0.4); hold on;
                h = fill([(1:10:trial_timepoints2)/100, trial_timepoints2/100],[parameters.pdgLUT{parameters.paradigms(k)}(1:10:trial_timepoints2,3); 0],'b');
                set(h, 'EdgeColor','none', 'FaceAlpha', 0.4);
                axis([0 trial_timepoints/7.5 0 3])
                plot((1:trial_timepoints)/7.5+50/7.5,smooth(parameters.plane{plane_nb}.timetraces{k}(1:trial_timepoints,neuron_nb),parameters.smoothing)/70,'k');
                hold off;

                % plot ROIs / neuron anatomy
                windowsize = 30;
                ROIX = squeeze(parameters.plane{plane_nb}.ROI_map(1,:,:));
                [x,y] = find(ROIX == neuron_nb); x = round(mean(x)); y = round(mean(y));
                xxx = max(1,x-windowsize):min(512,x+windowsize); yyy = max(1,y-windowsize):min(477,y+windowsize);
                ax(k+2*parameters.nb_trials) = subplot(nb_y,nb_x,k+ (row_pos)*nb_x+nb_x); imagesc(parameters.plane{plane_nb}.anatomy(xxx,yyy,k),[-30 70]); colormap(gray)
                title(parameters.FileNames(k,:),'Interpreter','None');
            end
            figure(918),
    end
end 