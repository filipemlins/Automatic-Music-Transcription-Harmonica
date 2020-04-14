function [plca_bin, plca_bin_s] = pitch_refinement(int_cqt, plca_out, thr_detection)

w = 3;
ref_note = 60;
bins_div = 5;

ss = size(plca_out);
plca_bin = zeros(ss);
plca_bin_s = zeros([round(ss(1)/5)+60, ss(2)]);
for i=1:ss(2);
    spec = int_cqt(:,i);    
    pitches = plca_out(:,i);    
    idx = find(pitches>thr_detection);    
    for j=1:length(idx);
            idx(j) = adjust_freq_to_peak(spec, idx(j), w); %% update pitch detections
    end
    plca_bin(idx,i)=1;
    %plca_bin_s(round((idx-60)/5+60 - 0.5),i)=1; %%TODO why -0.5?? adjusting the gt discretization
    plca_bin_s(round((idx-ref_note)/bins_div+ref_note),i)=1; %%TODO why -0.5?? adjusting the gt discretization
end;
