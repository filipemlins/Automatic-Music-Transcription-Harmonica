function [gt, gt2, gts] = create_gt_matrix(notes, plca_out, Fs, audio_length)

bins_div=5;
ref_note= 60; % not sure if the ref note is C4=60 or it is the lowest note in the dictionary //TODO double check this
s = size(plca_out);
ftime = (audio_length/Fs)/s(2);
gt = zeros(size(plca_out));
gt2 = zeros(size(plca_out));
gts = zeros(size(plca_out));

midi_pitch = round(freq2midi(notes(:,3)));
notes(:,1) = round(notes(:,1)./ftime);
notes(:,2) = round(notes(:,2)./ftime);
notes(:,3) = (freq2midi(notes(:,3))-ref_note)*bins_div+ref_note;

for i=1:size(notes,1)
    p = round(notes(i,3));
    if (p<1); continue; end;
    p2 = midi_pitch(i);
    gt(p,notes(i,1):notes(i,2))=1;
    gt2(p-1:p+1,notes(i,1):notes(i,2))=1;
    
    gts(p2,notes(i,1):notes(i,2))=1;    
end



