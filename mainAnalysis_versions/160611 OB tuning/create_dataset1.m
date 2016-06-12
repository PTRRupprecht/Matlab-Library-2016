% k = 1;

odors = {'Phe', 'Trp', 'His', 'Ala', 'Lys', 'Arg', 'Food'};
% sequence of odors:
odorIX = [4 4 6 6 1 1];

% when was the on-/offset of the odor valve
clear onset0 offset0
onset0(1:length(odorIX)) = 20;
offset0(1:length(odorIX)) = 55;
onset0([]) = NaN;
offset0([]) = NaN;

% additional delay due to flow in tube and discarding first frames
delay4 = 50/7.5-1.5;
delay8 = 50/7.5*2-1.5; % for 8 planes (lower framerate)

delay = delay4;

% write metadata of dataset, together with data from dataset 'plane'
% final structure containing everything is 'dataset1'
for jj = 1:numel(plane{1}.timetraces)
    dataset1{k}.dataset = '160608_OB2_pos3';
    dataset1{k}.odor = odors(odorIX(jj));
    dataset1{k}.comment = '';
    dataset1{k}.meta = plane{1}.meta;
    for ii = 1:numel(plane)
        dataset1{k}.ROI_map{ii} = squeeze(plane{ii}.ROI_map(jj,:,:));
        dataset1{k}.anatomy{ii} = squeeze(plane{ii}.anatomy(:,:,jj));
        dataset1{k}.timetraces_raw{ii} = plane{ii}.timetraces_raw{jj};
    end
    dataset1{k}.onset = onset0(jj) - delay;
    dataset1{k}.offset = offset0(jj) - delay;
    k = numel(dataset1) + 1; % increment counter for dataset1
end


