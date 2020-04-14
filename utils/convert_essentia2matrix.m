function m = convert_essentia2matrix(gt_mat, audio_length_sec, frame_list)

frame_tstep = audio_length_sec/size(gt_mat,2);
frame_idx = round(frame_list(:,1)/frame_tstep);

%frame_idx,
m = zeros(size(gt_mat));
p = round(freq2midi(frame_list(:,2)));
p(p>size(gt_mat,1),2) = 0;

frame_idx(frame_idx>size(gt_mat,2),:)=0;
list_p = [frame_idx, p];



for i=1:size(list_p,1)
    if list_p(i,1)>0 && list_p(i,2)>0 
        m(list_p(i,2), list_p(i,1)) = 1;    
    end
end







