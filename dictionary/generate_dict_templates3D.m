
%% this script extracts templates for the building task of the spectral dictionary
% please, observe that the audio file contains a sequence of pitches. Also, it is needed the
% respective pitch tracking from the YIN/pYIN plugin. We have used SonicVisualizer to perform the F0 extraction

%minFreq = 130.8128;  % midi = 48
minFreq = 27.5;
%bins_div = 1;  % bins per semitone
%bins_div = 3;  % bins per semitone
bins_div = 5;  % bins per semitone

wfile_l{1}= 'datasets/dictionary/harmonica/dict_harmonica_2Bb.wav';
wfile_l{2}= 'datasets/dictionary/harmonica/dict_harmonica_4g.wav';
wfile_l{3}= 'datasets/dictionary/harmonica/dict_harmonica_1c.wav';
wfile_l{4}= 'datasets/dictionary/harmonica/dict_harmonica_2c.wav';
wfile_l{5}= 'datasets/dictionary/harmonica/dict_harmonica_3c.wav';


K_KLUSTERS = 1;

out_dict_file = ['dict_harm',num2str(K_KLUSTERS),'f5d5_r.mat'];
%out_dict_file = 'dict_harm3.mat';
%out_dict_file = 'dict_harm5c.mat';
%out_dict_file = 'dict_harm5b.mat';
%out_dict_file = 'dict_harm20c.mat';
%out_dict_file = 'dict_harm10.mat';

dict5d = [];
for vd=1:length(wfile_l)

    wfile = wfile_l{vd};
if 1
    addpath('utils');
    addpath('cqt');
    
    %% load input wav file
    [y,Fs] = audioread(wfile);
    y=y(:,1);
    %% load ground truth (f0 seq)
    f0 = dlmread([wfile(1:end-4), '_f0.txt']);
    
    %% compute CQT
    %     maxFreq = Fs/2;
    %     bins = 12*bins_div;  % bins per octave
    %     sparKernel= sparseKernel(minFreq, maxFreq, bins, Fs);
    %     intCQT = schramm_cqt(y,sparKernel);
    %     intCQT(end:end-1,:)=0;
    %     intCQT = circshift(intCQT,[-2,0]); % TODO why? %
    % normalise by signal amplitude (it aims to remove background noise)
    %    intCQT=abs(intCQT);
    fmin = 27.5;
    [intCQT, Xcq] = computeVQT(wfile, 30);
    %[intCQT, Xcq] = computeCQT(wfile);
    intCQT(1:150,:)=0;
    
    
    sumx = sum(intCQT);
    %sumx(sumx<3)=0;
    sumx(sumx<0.005)=0;
    
    %     ss = diag(1./(sumx+eps));
    %     intCQT = intCQT*ss;
    %
    %     ss = diag((sumx>0.01));
    %     intCQT2 = intCQT*ss;
    %     figure;
    %     subplot(1,2,1); imagesc(intCQT); axis xy;
    %     subplot(1,2,2); imagesc(intCQT2); axis xy;
    %
    %     intCQT = intCQT2;
    %     sumx = sum(intCQT);
    
    %% remove negative freqs (from yin)
    f0(f0(:,2)<0,:) = [];
    
    %% adjust f0 time scale
    %tstep = 512/Fs;
    tstep = ((length(y)/Fs)/size(intCQT,2));
    stime = f0(:,1)./tstep;  %% convert pyin time frame to CQT time tsteps
    stime = round(stime);
    
    %% remove pitch detections on very low amplitudes (silence intervals)
    amp = find(sumx < 0.001);
    amp = amp(amp<=max(stime));
    for i=1:length(amp)
        f0(stime==amp(i),2) = -1;
    end
    f0(f0(:,2)<0,:) = [];
    
    %% adjust f0 time scale again since we have removed some frames
    stime = f0(:,1)./tstep;  %% convert pyin time frame to CQT time tsteps
    stime = round(stime);
    
end

%% mapping of input frequency to midi scale (and VQT bins)
% ref_note = 60; %% midi pitch reference (C4)
% lower_note = freq2midi(minFreq);
% m = (freq2midi(f0(:,2)) - ref_note) * bins_div + (ref_note-lower_note)*bins_div - 1;
lower_note = 21;
msemi = round(freq2midi(f0(:,2)));
m = (freq2midi(f0(:,2))-lower_note)*5 +1;

fig1 = figure;
imagesc(intCQT); axis xy; hold on;
plot(stime, m, 'r.');

%% detect onsets and offsets (based on silence gaps!)
th=20; %% TODO -- I am splitting the audio signal based on the onset gap from the f0 estimates using pyin
p = find(abs(diff(stime))>=th); %% get the index of each f0 gap (onsets)
pEnd = [p; length(diff(stime))];
pIni = [1; p+1];

%% alloc structures
dict = [];
pvec = [];
%ener = [];

%% plotting
%s = ones([1 length(p)])*size(sparKernel,2);
s = ones([1 length(p)])*size(intCQT,1);
stem(stime(p),s,'m');
fig2 = figure;
pause(.1),

%% scan pitch (f0) detections and extract the spectral template for each of them
% (also performs the pitch shift over 20 cent resolution)
%iniP = 1;
prevPitch = [];
for i=1:length(pIni);
    
    figure(fig1);
    stem(stime(pIni(i)), size(intCQT,1), 'g');
    stem(stime(pEnd(i)), size(intCQT,1), 'y');
    
    pitches = m(pIni(i):pEnd(i));
    pitches_semi = msemi(pIni(i):pEnd(i));
    %if (std(pitches) > bins_div) continue; end; %% avoid spurious yin detections
    if (length(pitches) <25) continue; end; %% avoid spurious yin detections
    
    pitch = round(median(pitches)); % remove outliers
    pitches(abs(pitches-pitch)>bins_div)=[];
    
    pitch = round(median(pitches)); % /TODO make it better
    pitch_semi = round(median(pitches_semi));
    
    %pitch = round(median(pitches(10:end-10))); % /TODO make it better
    time_pitches = round(stime(pIni(i):pEnd(i)));
    
    
    
        if 1
            %% fill the gaps between semitones
            if ~isempty(prevPitch)
                c=bins_div;
                shifted_tpl = dict(prevPitch,:);
                for j=prevPitch+1:pitch_semi-1
                    dict(j,:) = circshift(shifted_tpl,[0 c]);
                    pvec(j) = j ; % please, double check this freq offset
                    c=c+bins_div;
                end
            end
        end
    
    
    
    
    
    seg = intCQT(:,time_pitches);
    
    %figure; imagesc(seg);
    sss = [];
    for f=1:size(seg,2)
        [~, shift_step] = adjust_freq_to_peak(seg(:,f), pitch, 4);
        sss = [sss, shift_step];
        %seg(:,f) = circshift(seg(:,f), -shift_step);
    end
    %sss,
    seg(:,sss~=0)=[];
    
    if size(seg,2)<=5
        continue;
    end
    %figure; imagesc(seg);
    IX = kmeans(seg', K_KLUSTERS); %% five clusters
    
    for k=1:K_KLUSTERS
        tpl_k = median(seg(:,IX==k),2);
        %     tpl2 = median(seg(:,IX==2),2);
        %     tpl3 = median(seg(:,IX==3),2);
        %     tpl4 = median(seg(:,IX==4),2);
        %     tpl5 = median(seg(:,IX==5),2);
        
        %tpl = median(seg,2);
        dict(pitch_semi,:,k) = tpl_k;
        %     dict(pitch_semi,:,2) = tpl2;
        %     dict(pitch_semi,:,3) = tpl3;
        %     dict(pitch_semi,:,4) = tpl4;
        %     dict(pitch_semi,:,5) = tpl5;
    end
    
    pvec(pitch_semi) = pitch_semi;
    prevPitch = pitch_semi;
    
    
    figure(fig2);
    subplot(2,1,1); imagesc(dict(pitch_semi,:,1)); axis xy;
    subplot(2,1,2); plot(dict(pitch_semi,:,k) );
    title(sprintf('%2.2f',pitch_semi));
    pause(.1),
end

% making the 3d shifted tensor

dict= permute(dict, [2 1 3]);

hc = createHarmonicComb();
for i=1:K_KLUSTERS
    dict(:,:,i) = dict(:,:,i).*hc(1:545,1:size(dict,2));
end

% normalise templates
for i=1:size(dict,2);
    for j=1:size(dict,3);
        dict(:,i,j) = dict(:,i,j)./(sum(dict(:,i,j))+eps);
    end
end

shifted_tpl = dict(:,60);
c=bins_div;
for i=59:-1:54    
  dict(:,i) = circshift(shifted_tpl,[-c 0]);
  c=c+bins_div;
end

shifted_tpl = dict(:,91);
c=bins_div;
for i=92:96    
  dict(:,i) = circshift(shifted_tpl,[c 0]);
  c=c+bins_div;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
f1 = circshift(dict, [-2, 0, 0]);
f2 = circshift(dict, [-1, 0, 0]);
f3 = dict;
f4 = circshift(dict, [1, 0, 0]);
f5 = circshift(dict, [2, 0, 0]);

dict = [];
dict(:,:,:,1) = f1;
dict(:,:,:,2) = f2;
dict(:,:,:,3) = f3;
dict(:,:,:,4) = f4;
dict(:,:,:,5) = f5;


dict5d(:,:,:,:,vd) = permute(dict, [1 2 4 3]);
sd5d = dict5d(:,54:end,:,:,:);

end
sd5d(:,1,:,:) = 0; % first pitch bin is used to hold all spurious activations when audio signal is in silence

save(out_dict_file, 'dict', 'dict5d', 'sd5d', 'bins_div', 'minFreq', 'lower_note',  'pvec', 'intCQT');

%figure; imagesc(dict); axis xy;


