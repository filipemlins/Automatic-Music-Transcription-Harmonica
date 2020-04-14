%% Demonstration
%Illustration of the use of programs to compute the
%forward and inverse constant-Q transform
%
%Christian Schörkhuber, Anssi Klapuri 2010-06

%% init values for CQT
fs = 44100;
bins_per_octave = 24;
fmax = fs/3;     %center frequency of the highest frequency bin 
fmin = fmax/512; %lower boundary for CQT (lowest frequency bin will be immediately above this): fmax/<power of two> 

%% generate/read input signal%
x = randn(30*fs,1);
% Drop frequencies outside [fmin fmax] to allow calculating 
% the SNR after inverse transform
x = [zeros(500,1); x; zeros(500,1)];
w1 = 2*(fmin/(fs/2)); w2 = 0.8*(fmax/(fs/2));
[B,A] = butter(6,[w1 w2]); x = filtfilt(B,A,x); 

%% CQT
Xcqt = cqt(x,fmin,fmax,bins_per_octave,fs);

%***computing cqt with optional input parameters***********
% Xcqt = cqt(x,fmin,fmax,bins_per_octave,fs,'q',1,'atomHopFactor',0.25,'thresh',0.0005,'win','sqrt_blackmanharris');

%***computing rasterized complex coefficients***************
% Xcqt = cqtPerfectRast(x,fmin,fmax,bins_per_octave,fs,'q',1,'atomHopFactor',0.25,'thresh',0.0005,'win','sqrt_blackmanharris');

%***precomputing filter and kernel**************************
% [B A] = butter(6,0.5,'low'); %design f_nyquist/2-lowpass filter
% K = genCQTkernel(fmax,bins_per_octave,fs,'atomHopFactor',atomHopFactor);
% Xcqt = cqt(x,fmin,fmax,bins_per_octave,fs,'kernel',K,'coeffB',B,'coeffA',A);

%***precomputing filter and kernel using cqtPerfectRast*****
% [B A] = butter(6,0.5,'low'); %design f_nyquist/2-lowpass filter
% K = genCQTkernel(fmax,bins_per_octave,fs,'atomHopFactor',atomHopFactor,'perfRast',1);
% Xcqt = cqtPerfectRast(x,fmin,fmax,bins_per_octave,fs,'kernel',K,'coeffB',B,'coeffA',A);

%% inverse CQT
y = icqt(Xcqt);
SNR = 10*log10(sum(abs(x).^2)/sum(abs(x-y).^2)); %calculate the SNR

%% plot CQT
 plotCQT(Xcqt,fs,0.6,'surf');

%% get CQT (interpolated)
% intCQT = getCQT(Xcqt,'all','all');

