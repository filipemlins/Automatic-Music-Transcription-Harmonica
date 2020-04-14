function notes = plot_gt_from_file(gt_txt_file, pitch_scale)

if nargin < 2
    pitch_scale = 0;
else
    switch lower(pitch_scale)
        case {'midi', 'log'}
            pitch_scale = 1;
        case 'hz'
            pitch_scale = 0;
        otherwise
            pitch_scale = 0;
    end
end

notes = dlmread(gt_txt_file);
s = size(notes);

if notes(1,1)-notes(1,2)>0
    notes(:,2) = notes(:,1)+notes(:,2);
end

if pitch_scale   %% == midi
    notes(:,3) = midi2freq(notes(:,3));
end

hold on;
for i=1:s(1);
    line([notes(i,1), notes(i,2)], freq2midi([notes(i,3),notes(i,3)]), 'LineWidth', 2, 'Color', 'm');
end