

function getF0(num_frames, y,Fs)

f0 = yin(y,Fs);

f = f0.f0*440;
f(isnan(f))=0;
f(f<0)=0;

f = freq2midi(f);
f(f<0)=0;
f(isnan(f))=0;
f = f*5+21;

frame_idx = round( (1:length(f))/length(f) * num_frames);

figure; imagesc(intCQT); axis xy; hold on
plot(frame_idx, f);


%frame_idx,
m = zeros(size(gt_mat));
p = round(freq2midi(frame_list(:,2)));
p(p>size(gt_mat,1),2) = 0;

frame_idx(frame_idx>size(gt_mat,2),:)=0;
list_p = [frame_idx, p];