function [intVQT,Xcq] = computeVQT(filename)
% Settings for computing VQT for music signals (by E. Benetos)
%addpath('VQT');

% PARAMETERS
fs = 44100;
fmin = 27.5;
B = 60;
gamma = 30; 
fmax = fs/2;


% Load .wav file
[x,fs] = audioread(filename);
if (size(x,2) == 2) y = mean(x')'; else y=x; end;
if (fs ~= 44100) y = resample(y,44100,fs); end;
fs = 44100;


% Compute VQT
Xcq = vqt(y, B, fs, fmin, fmax, 'rasterize', 'full','gamma', gamma);
intVQT = Xcq.c;
intVQT = abs(intVQT(1:545,:));