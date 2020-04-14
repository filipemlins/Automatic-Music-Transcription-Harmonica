
%% this script extracts templates for the building task of the spectral dictionary
% please, observe that the audio file contains a sequence of pitches. Also, it is needed the
% respective pitch tracking from the YIN/pYIN plugin. We have used SonicVisualizer to perform the F0 extraction

minFreq = 130.8128;  % midi = 48    
%bins_div = 1;  % bins per semitone
%bins_div = 3;  % bins per semitone
bins_div = 5;  % bins per semitone
wfile= 'datasets/dictionary/harmonica/dict_harmonica_1c.wav';
%wfile= 'datasets/dictionary/harmonica/dict_harmonica_2c.wav';
%wfile= 'datasets/dictionary/harmonica/dict_harmonica_3c.wav';

%out_dict_file = 'dict_harm1.mat';
%out_dict_file = 'dict_harm3.mat';
out_dict_file = 'dict_harm5b.mat';
%out_dict_file = 'dict_harm10.mat';


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
% normalise by signal amplitude (it aims to remove background noise)
%    intCQT=abs(intCQT);
        fmin = 27.5;
        [intCQT, Xcq] = computeVQT(input_wav, 30);

    sumx = sum(intCQT);
    ss = diag(1./(sumx+eps));
    intCQT = intCQT*ss;
    
    ss = diag((sumx>0.01));
    intCQT2 = intCQT*ss;
    figure;
    subplot(1,2,1); imagesc(intCQT); axis xy;
    subplot(1,2,2); imagesc(intCQT2); axis xy;
    
    intCQT = intCQT2;
    sumx = sum(intCQT);
    
    %% remove negative freqs (from yin)
    f0(f0(:,2)<0,:) = [];
    
    %% adjust f0 time scale
    step = 512/Fs;
    stime = f0(:,1)./step;  %% convert pyin time frame to CQT time steps
    stime = round(stime);
    
    %% remove pitch detections on very low amplitudes (silence intervals)
    amp = find(sumx < 0.1);
    amp = amp(amp<=max(stime));
    for i=1:length(amp)
        f0(stime==amp(i),2) = -1;
    end
    f0(f0(:,2)<0,:) = [];
    
    %% adjust f0 time scale again since we have removed some frames
    stime = f0(:,1)./step;  %% convert pyin time frame to CQT time steps
    stime = round(stime);
    
end

%% mapping of input frequency to midi scale (and VQT bins)
ref_note = 60; %% midi pitch reference (C4)
lower_note = freq2midi(minFreq);
m = (freq2midi(f0(:,2)) - ref_note) * bins_div + (ref_note-lower_note)*bins_div + 1;

fig1 = figure;
imagesc(intCQT2); axis xy; hold on;
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
s = ones([1 length(p)])*size(sparKernel,2);
stem(stime(p),s,'m');
fig2 = figure;
pause(.1),

%% scan pitch (f0) detections and extract the spectral template for each of them
% (also performs the pitch shift over 20 cent resolution)
%iniP = 1;
prevPitch = [];
for i=1:length(pIni);
    
    figure(fig1);
    stem(stime(pIni(i)), size(sparKernel,2), 'g');
    stem(stime(pEnd(i)), size(sparKernel,2), 'y');
    
    pitches = m(pIni(i):pEnd(i));
    if (std(pitches) > bins_div) continue; end; %% avoid spurious yin detections
    
    
    pitch = round(median(pitches)); % /TODO make it better
    %pitch = round(median(pitches(10:end-10))); % /TODO make it better
    time_pitches = round(stime(pIni(i):pEnd(i)));
    
    
    
    if 1
        %% fill the gaps between semitones
        if ~isempty(prevPitch)
            c=1;
            shifted_tpl = dict(prevPitch,:);
            for j=prevPitch+1:pitch-1
                dict(j,:) = circshift(shifted_tpl,[0 c]);
                pvec(j) = j ; % please, double check this freq offset
                c=c+1;
            end
        end
    end
    
    seg = intCQT(:,time_pitches);
    tpl = median(seg,2);
    dict(pitch,:) = tpl;
    pvec(pitch) = pitch;
    prevPitch = pitch;
    
    
    figure(fig2);
    subplot(2,1,1); imagesc(seg); axis xy;
    subplot(2,1,2); plot(tpl);
    title(sprintf('%2.2f',pitch));
    pause(.1),
end

%% TODO
% harmonicComb filter

dict= permute(dict, [2 1 3]);
sd = dict(:,round((ref_note-lower_note)*bins_div)+1:end);

save(out_dict_file, 'dict', 'sd', 'bins_div', 'minFreq', 'lower_note', 'ref_note', 'pvec');

figure; imagesc(dict); axis xy;


