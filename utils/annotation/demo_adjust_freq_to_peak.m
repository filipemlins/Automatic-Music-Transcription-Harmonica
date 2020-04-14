minFreq = 130.8128;  % midi = 48
%bins_div = 1;  % bins per semitone
%bins_div = 3;  % bins per semitone
bins_div = 5;  % bins per semitone
input_wav = 'datasets/teste01/heyjudeG01.wav';
input_gt = 'datasets/teste01/heyjudeG01.txt';

%load  dict_harm1.mat;
%load  dict_harm3.mat;
load  dict_harm5.mat;


addpath('utils');
addpath('cqt');

%% load input wav file
[y,Fs] = audioread(input_wav);
y=y(:,1);

%% compute CQT
maxFreq = Fs/2;
bins = 12*bins_div;  % bins per octave
sparKernel= sparseKernel(minFreq, maxFreq, bins, Fs);
intCQT = schramm_cqt(y,sparKernel);


figure; imagesc(intCQT); axis xy;
time_frame = 200;
freq_pos = 207; %% just a guess (you must replace this by your iterative user pixel selection)
w = 10; 
%% look at the closest peak (at time frame 200)
new_freq_pos = adjust_freq_to_peak(intCQT(:,time_frame), freq_pos, w);

new_freq_pos,