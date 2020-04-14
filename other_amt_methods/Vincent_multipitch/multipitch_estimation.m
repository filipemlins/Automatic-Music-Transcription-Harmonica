function multipitch_estimation(wavfile,f0file)

% MULTIPITCH_ESTIMATION Frame-by-frame multiple pitch estimation on the
% MIDI scale using NMF under harmonicity and spectral smoothness constraints
%
% multipitch_estimation(wavfile,f0file)
%
% Inputs:
% wavfile: input audio filename
% f0file: output transcription filename (MIREX format with discrete pitches
% in Hz every 10 ms)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2009 Emmanuel Vincent
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
% If you find it useful, please cite the following reference:
% Emmanuel Vincent, Nancy Bertin and Roland Badeau, "Adaptive harmonic
% spectral decomposition for multiple pitch estimation," IEEE Trans. on
% Audio, Speech and Language Processing, to appear.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Defining NMF parameters
nbfreq=250;
nbcomp=1;
maxclus=6;
totwidth=22;
clusspace=totwidth/maxclus;
cluswidth=clusspace*2;
beta=.5;
thresh=27;

% Defining minimum note duration
mindur=eps;

% Reading WAV file
[x,fs]=wavread(wavfile);
x=resample(x,22050,fs).';

% Applying NMF
[U,A]=nmf_harmclus_erbt(x,22050,nbfreq,nbcomp,maxclus,clusspace,cluswidth,beta);

% Estimating onset/offset/f0
testnotes=harm2midi_erbt(U,A,22050,thresh,mindur);
testnotes(:,3)=440*2.^((testnotes(:,3)-69)/12);

% Transcribing into frame-by-frame format
N=ceil(length(x)/220.5);
test=zeros(N,1);
for note=1:size(testnotes,1),
    b=max(1,floor(testnotes(note,5)*100)+1);
    e=min(N,floor(testnotes(note,6)*100)+1);
    p=1;
    while any(test(b:e,p)),
        p=p+1;
        if p > size(test,2),
            test=[test zeros(N,1)];
        end
    end
    test(b:e,p)=testnotes(note,3);
end
test=test(1:N,:);

% Writing into a file
fid=fopen(f0file,'w');
for n=1:N-1,
    fprintf(fid,'%.2f\t',(n-1)/100);
    testloc=test(n,test(n,:)~=0);
    nbloc=length(testloc);
    for p=1:nbloc-1,
        fprintf(fid,'%.2f\t',testloc(p));
    end
    if nbloc,
        fprintf(fid,'%.2f\n',testloc(nbloc));
    else
        fprintf(fid,'\n');
    end
end
fprintf(fid,'%.2f\t',(N-1)/100);
testloc=test(N,test(N,:)~=0);
nbloc=length(testloc);
for p=1:nbloc-1,
    fprintf(fid,'%.2f\t',testloc(p));
end
if nbloc,
    fprintf(fid,'%.2f',testloc(nbloc));
end
fclose(fid);

return;



function [U,A,dist,B]=nmf_harmclus_erbt(x,fs,F,nbcomp,maxclus,clusspace,cluswidth,beta)

% NMF_HARMCLUS_ERBT Harmonic NMF of a signal from a ERB transform with
% basis spectra representing partial clusters and fundamental frequencies
% on the MIDI scale tuned at 440 Hz
%
% Initialization with a slope of 6 dB/octave for the first component, 12
% dB/oct for the second, etc
%
% [U,A,dist,B]=nmf_harmclus_erbt(x,fs,F,nbcomp,maxclus,clusspace,cluswidth,beta)
%
% Inputs:
% x: 1 x T vector containing a single-channel signal
% fs: sampling frequency in Hz
% F: number of frequency bins
% nbcomp: number of spectral envelope components per note
% maxclus: maximal number of partial clusters per note
% clusspace: spacing between successive partial clusters in ERB
% cluswidth: bandwidth of each partial cluster in ERB
% beta: distortion measure (0-> IS, 1-> KL, 2-> EUC)
%
% Output:
% U: F x (nbnotes x nbcomp) matrix containing NMF basis vectors
% A: (nbnotes x nbcomp) x N matrix containing NMF time weights
% dist: achieved distortion measure
% B: nbclus x nbcomp matrix containing spectral envelopes

%%% Errors and warnings %%%
if nargin<8, error('Not enough input arguments.'); end
[I,T]=size(x);
if I>T, error('The input signal must contain more time samples than channels.'); end
if fs>25000, error('The sampling frequency must be smaller than 25 kHz.'); end
wlen=2^nextpow2(.02*fs);    %20 ms window length
N=ceil(T/wlen);

%%% Computing ERBT coefficients and frequency scale %%%
X=zeros(F,N,I);
for i=1:I,
    [X(:,:,i),f]=erbtm(x(i,:),fs,F,wlen);
end
X=(sum(X.^2,3)+1e-18).^.5;
fmin=f(1); fmax=f(F);
emin=9.26*log(.00437*fmin+1); emax=9.26*log(.00437*fmax+1);
e=(0:F-1)*(emax-emin)/(F-1)+emin;
a=.5*(F-1)/(emax-emin)*9.26*.00437*fs*exp(-e/9.26)-.5;
alen=2*round(a)+1;
f=f/fs;

%%% Defining partial spectra %%%
firstnote=21; lastnote=108; nbnotes=lastnote-firstnote+1; pitch=firstnote:lastnote;
f0=2.^((pitch-69)/12)*440/fs;
nharm=floor(.5./f0);
ppos=[0,cumsum(nharm)];
nbpart=ppos(end);
Z=zeros(F,nbpart);   %partial spectra
partfreq=zeros(1,nbpart);
for n=1:nbnotes,
    partfreq(ppos(n)+1:ppos(n)+nharm(n))=f0(n)*(1:nharm(n));
end
for c=1:F,
    Z(c,:)=abs(sinc((f(c)-partfreq)*alen(c))+.5*sinc((f(c)-partfreq)*alen(c)+1)+.5*sinc((f(c)-partfreq)*alen(c)-1));
end

%%% Defining cluster spectra and initializing note spectra %%%
clusnum=zeros(1,nbpart);    % cluster index for each partial
for n=1:nbnotes,
    clusnum(ppos(n)+1:ppos(n)+nharm(n))=9.26*(log(.00437*partfreq(ppos(n)+1:ppos(n)+nharm(n))*fs+1)-log(.00437*partfreq(ppos(n)+1)*fs+1))/clusspace;
end
nclus=min(maxclus,round(clusnum(ppos(2:end)))+1);
cpos=[0,cumsum(nclus)];
nbclus=cpos(end);
clusfreq=zeros(1,nbclus);   % center frequency of each cluster
for n=1:nbnotes,
    clusfreq(cpos(n)+1:cpos(n)+nclus(n))=((.00437*partfreq(ppos(n)+1)*fs+1)*exp(clusspace/9.26*(0:nclus(n)-1))-1)/(.00437*fs);
end
V=zeros(F,nbclus);   %cluster spectra
B=zeros(nbclus,nbcomp);  %cluster weights
U=zeros(F,nbnotes*nbcomp);  %note spectra
for n=1:nbnotes,
    weights=zeros(nharm(n),nclus(n));
    for c=1:nclus(n),
        r=(clusnum(ppos(n)+1:ppos(n)+nharm(n))-c+1)*clusspace/cluswidth;
        order=4;
        k=sqrt(pi)*gamma(order-.5)/gamma(order);
        weights(:,c)=(1+(k*r).^2).^-order;
    end
    V(:,cpos(n)+1:cpos(n)+nclus(n))=Z(:,ppos(n)+1:ppos(n)+nharm(n))*weights;
    B(cpos(n)+1:cpos(n)+nclus(n),:)=10.^(-6/20*log2(clusfreq(cpos(n)+1:cpos(n)+nclus(n))/clusfreq(cpos(n)+1)).'*(1:nbcomp));
    U(:,(0:nbcomp-1)*nbnotes+n)=V(:,cpos(n)+1:cpos(n)+nclus(n))*B(cpos(n)+1:cpos(n)+nclus(n),:);
end

%%% Performing NMF updates %%%
A=ones(nbnotes*nbcomp,N);
Y=U*A;
gconverged=0; dist=inf;
while ~gconverged,
    gprevdist=dist;
    % Updating note weights
    lconverged=0;
    while ~lconverged,
        lprevdist=dist;
        A=A.*(U.'*(X.*Y.^(beta-2)))./(U.'*(Y.^(beta-1)));
        Y=U*A;
        switch beta
            case 0,
                dist=sum(sum(X./Y-log(X./Y)-1));
            case 1,
                dist=sum(sum(X.*log(X./Y)+Y-X));
            case 2,
                dist=.5*sum(sum((X-Y).^2));
            otherwise
                dist=sum(sum(X.^beta+(beta-1)*Y.^beta-beta*X.*Y.^(beta-1)))/(beta*(beta-1));
        end
        lconverged=(10*log10(lprevdist/dist) < 5e-3);
    end
    % Updating note spectra
    lconverged=0;
    while ~lconverged,
        lprevdist=dist;
        XYA=(X.*Y.^(beta-2))*A.';
        YA=(Y.^(beta-1))*A.';
        for n=1:nbnotes,
            B(cpos(n)+1:cpos(n)+nclus(n),:)=B(cpos(n)+1:cpos(n)+nclus(n),:).*(V(:,cpos(n)+1:cpos(n)+nclus(n)).'*XYA(:,(0:nbcomp-1)*nbnotes+n))./(V(:,cpos(n)+1:cpos(n)+nclus(n)).'*YA(:,(0:nbcomp-1)*nbnotes+n,:)+realmin);
            U(:,(0:nbcomp-1)*nbnotes+n)=V(:,cpos(n)+1:cpos(n)+nclus(n))*B(cpos(n)+1:cpos(n)+nclus(n),:);
        end
        Y=U*A;
        switch beta
            case 0,
                dist=sum(sum(X./Y-log(X./Y)-1));
            case 1,
                dist=sum(sum(X.*log(X./Y)+Y-X));
            case 2,
                dist=.5*sum(sum((X-Y).^2));
            otherwise
                dist=sum(sum(X.^beta+(beta-1)*Y.^beta-beta*X.*Y.^(beta-1)))/(beta*(beta-1));
        end
        lconverged=(10*log10(lprevdist/dist) < 5e-3);
    end
    % Convergence test
    gconverged=(10*log10(gprevdist/dist) < 1e-2);
end

%%% Energy normalization %%%
A=A.*(sum(U.^2).^.5.'*ones(1,N));
U=U./(ones(F,1)*sum(U.^2).^.5);

return;



function [X,f]=erbtm(x,fs,F,wlen)

% ERBTM Magnitude ERB Transform using a Hann window.
%
% [X,f]=erbtp(x,fs,F,wlen)
%
% Inputs:
% x: 1 x T vector containing a single-channel signal
% fs: sampling frequency in Hz
% F: number of frequency bins (the ratio between the bandwidth of each bin
% and the frequency difference between successive bins is constant)
% wlen: number of samples per frame (must be a multiple of the largest
% downsampling factor, typically a large power of 2)
%
% Output:
% X: F x N matrix containing the time-frequency magnitude (amplitude) coefficients
% f: F x 1 vector containing the center frequency of each frequency bin

%%% Errors and warnings %%%
if nargin<4, error('Not enough input arguments.'); end
[I,T]=size(x);
if I>1, error('The input signal must be a row vector.'); end
N=ceil(T/wlen);

%%% Computing ERBT coefficients %%%
x=hilbert(x);
X=zeros(F,N);
% Determining minimum and maximum frequency
fmax=.5*fs; fmin=0;
for j=1:100,
    emin=9.26*log(.00437*fmin+1);
    emax=9.26*log(.00437*fmax+1);
    fmin=1.5*(emax-emin)/(F-1)/9.26/.00437*exp(emin/9.26);
    fmax=.5*fs-1.5*(emax-emin)/(F-1)/9.26/.00437*exp(emax/9.26);
    if (fmax < 0) || (fmin > .5*fs), error('The number of frequency bins is too small.'); end
end
% Determining frequency and window length scales
emin=9.26*log(.00437*fmin+1);
emax=9.26*log(.00437*fmax+1);
e=(0:F-1)*(emax-emin)/(F-1)+emin;
f=(exp(e/9.26)-1)/.00437;
a=.5*(F-1)/(emax-emin)*9.26*.00437*fs*exp(-e/9.26)-.5;
% Determining dyadic downsampling bins (for fast computation)
fup=f+1.5*fs./(2*a+1);
subs=-log(2*fup/fs)/log(2);
subs=2.^max(0,floor(min(log2(wlen),subs)));
if (wlen/subs(1) ~= floor(wlen/subs(1))), error(['The number of samples per frame must be a multiple of ' int2str(subs(1)) '.']); end
down=(subs~=[subs(2:end),1]);
for bin=F:-1:1,
    % Dyadic downsampling
    if down(bin),
        x=resample(x,1,2,50);
    end
    % Convolution with a modulated sine window
    hwlen=round(a(bin)/subs(bin));
    win=hanning(2*hwlen+1).'.*exp(2i*pi*f(bin)/fs*subs(bin)*(-hwlen:hwlen));
    band=[fftfilt(win,[x,zeros(1,2*hwlen)]) zeros(1,wlen/subs(bin))];
    % Square-root energy on short time frames
    band=band(hwlen+1:hwlen+N*wlen/subs(bin));
    X(bin,:)=sqrt(sum(reshape(abs(band).^2,wlen/subs(bin),N),1)/(hwlen+1)^2*subs(bin));
end

return;



function nmat=harm2midi_erbt(W,H,fs,thresh,mindur)

% HARM2MIDI_ERBT Processing of harmonic NMF outputs into a note matrix
%
% nmat=harm2midi_erbt(W,H,fs,thresh,mindur)
% 
% Inputs:
% W: F x (nbnotes x nbcomp) matrix containing NMF basis vectors
% H: (nbnotes x nbcomp) x N matrix containing NMF time weights
% fs : sampling rate
% thresh: relative onset amplitude threshold in dB
% mindur: minimum note duration in seconds
% 
% Output:
% nmat: midi-like matrix of the transcription result (Ken Schutte's MIDI
% Toolbox format)

n=size(H,2);
r=size(W,2);
wlen=2^nextpow2(.02*fs);

% Pitch definition
firstnote=21; lastnote=108; nbnotes=lastnote-firstnote+1;
nbcomp=r/nbnotes;
pitch=reshape((firstnote:lastnote).'*ones(1,nbcomp),1,nbnotes*nbcomp);

% Summation of components with the same pitch
k=1;
while k <= size(H,1)
    currpitch = pitch(k);
    indices = find(pitch(k+1:end)==currpitch);
    if ~isempty(indices)
        X=W(:,[k indices+k])*H([k indices+k],:);
        H(k,:) = sqrt(sum(X.^2));
        H(indices+k,:)=[];
        pitch(indices+k)=[];
    end
    k = k+1;
end
r = size(H,1);

% Minimal note magnitude and duration
thresh=max(max(H))*10^(-thresh/20);
mindur=ceil(mindur*fs/wlen);

% TH coefficients are 1 for a "note ON" event
TH = zeros(r,n);
TH(:,1) = all((H(:,1:mindur)>thresh),2);
for k=2:n-mindur+1
    TH(:,k) = (H(:,k-1)<=thresh) & all((H(:,k:k+mindur-1)>thresh),2);
end

% DH coefficients are 1 for a "note OFF" event
DH = zeros(r,n);
for k=mindur:n-1
    DH(:,k) = (H(:,k+1)<=thresh) & all((H(:,k-mindur+1:k)>thresh),2);
end
DH(:,n) = all((H(:,n-mindur+1:n)>thresh),2);

% initialize notematrix
pos_events = find(TH==1);
nb_events = length(pos_events);
nmat = zeros(nb_events,8);

% builds output matrix
for p=1:nb_events
    [atom,timeon]=ind2sub([r n],pos_events(p));
    if timeon > 1,
        timeon=timeon-(log(H(atom,timeon))-log(thresh))/(log(H(atom,timeon))-log(H(atom,timeon-1)));
    end
    timeoff=find(DH(atom,ceil(timeon):end)); timeoff=ceil(timeon)+timeoff(1)-1;
    if timeoff < n,
        timeoff=timeoff+(log(H(atom,timeoff))-log(thresh))/(log(H(atom,timeoff))-log(H(atom,timeoff+1)));
    end
    nmat(p,3)=pitch(atom);
    nmat(p,4)=100;
    nmat(p,5)=(timeon-.5)*wlen/fs;
    nmat(p,6)=(timeoff-.5)*wlen/fs;
end
[nmat(:,5),ind]=sort(nmat(:,5));
nmat(:,[1:4 6])=nmat(ind,[1:4 6]);
nmat(:,1)=1;
nmat(:,2)=0;
[time,ind]=sort([nmat(:,5);nmat(:,6)]);
for p=1:2*nb_events,
    nmat(6*nb_events+ind(p))=p;
end

% delete unrelevant events
nmat=nmat((nmat(:,3)>0),:);

return;