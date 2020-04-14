function Xcqt = cqtPerfectRast(x,fmin,fmax,bins,fs,varargin)
%Xcqt = cqtPerfectRast(x,fmin,fmax,bins,fs,varargin) calculates the
%constant-Q transform of the input signal x so that complex coefficients
%are computed at the time resolution of the hiqhest frequency for all octaves.
%This means that the cqt is actually performed several times for lower octaves 
%to provide also accurate phase information for the rasterized representation of
%the cqt coefficients (slower than cqt()).
%Note that if the kernel structure is computed prior to
%cqtPerfectRast() ('kernel' is an optional input parameter), the flag
%'perfRast' in the genCQTkernel() call has to be set to 1.
%
%INPUT:
%   fmin      ... lowest frequency of interest
%   fmax      ... highest frequency of interest
%   bins      ... frequency bins per octave
%   fs        ... sampling rate
%
%   optional input parameters (parameter name/value pairs):
%
%   'atomHopFactor' ... overlap of temporal atoms in percent. Default: 0.25.
%    
%   'q'             ... the maximum value for optimal reconstruction is q=1.
%                       For values smaller than 1 the bandwidths of the spectral
%                       atoms (filter) are increased retaining their center
%                       frequencies (frequency 'smearing', frequency domain redundancy 
%                       increases, time resolutin improves). Default: 1.
%   'thresh'        ... all values in the cqt kernel smaller than tresh are
%                       rounded to zero. A high value for thresh yields a
%                       very sparse kernel (fast) but introduces a bigger error. 
%                       The default value is chosen so that the error due to rounding is negligible.
%   'kernel'        ... if the cqt kernel structure has been precomputed
%                       (using function 'genCQTkernel'), the computation of the kernel
%                       will be by-passed below).
%   'win'           ... defines which window will be used for the CQT. Valid
%                       values are: 'blackman','hann' and 'blackmanharris'. To
%                       use the square root of each window use the prefix 'sqrt_'
%                      (i.e. 'sqrt_blackman'). Default: 'sqrt_blackmanharris'
%   'coeffB',
%   'coeffA'        ... Filter coefficients for the anti-aliasing filter, where
%                      'coeffB' is the numerator and 'coeffA' is the
%                       denominator (listed in descending powers of z). 
%                                                  
%OUTPUT:
%   Xcqt      ... struct that comprises various fields: 
%               spCQT: CQT coefficients in the form of a sparse matrix 
%                     (rasterized, not interpolated)
%               fKernel: spectral Kernel 
%               fmin: frequency of the lowest bin
%               fmax: frequency of the hiqhest bin
%               octaveNr: number of octaves processed
%               bins: number of bins per octave
%               intParams: structure containing additional parameters for the inverse transform   
%
%Christian Schörkhuber, Anssi Klapuri 2010-06
%2011-03: error removed in the calculation of fmin for output in Xcqt structure

%% input checking
if size(x,2) > 1 && size(x,1) > 1, error('cqt requires one-dimensional input!'); end;
if size(x,2) > 1, x = x(:); end; %column vector

%% input parameters
q = 1; %default value
atomHopFactor = 0.25; %default value
thresh = 0.0005; %default value
winFlag = 'sqrt_blackmanharris';

for ain = 1:2:length(varargin)
    if strcmp(varargin{ain},'q'), q = varargin{ain+1}; end;
    if strcmp(varargin{ain},'atomHopFactor'), atomHopFactor = varargin{ain+1}; end;
    if strcmp(varargin{ain},'thresh'), thresh = varargin{ain+1}; end;
    if strcmp(varargin{ain},'kernel'), cqtKernel = varargin{ain+1}; end;
    if strcmp(varargin{ain},'win'), winFlag = varargin{ain+1}; end;
    if strcmp(varargin{ain},'coeffB'), B = varargin{ain+1}; end;
    if strcmp(varargin{ain},'coeffA'), A = varargin{ain+1}; end;    
end

%% define
octaveNr = ceil(log2(fmax/fmin));
fmin = (fmax/2^octaveNr) * 2^(1/bins); %set fmin to actual value
xlen_init = length(x);

%% design lowpass filter
if ~exist('B','var') || ~exist('A','var')
    LPorder = 6; %order of the anti-aliasing filter
    cutoff = 0.5;
    [B A] = butter(LPorder,cutoff,'low'); %design f_nyquist/2-lowpass filter
end

%% design kernel for one octave
if exist('cqtKernel','var')
    if cqtKernel.perfRast == 0
        clear cqtKernel;
        warning('CQT:INPUT','Wrong kernel design! New kernel will be generated.');
    end
end
if ~exist('cqtKernel','var')      
    cqtKernel = genCQTkernel(fmax, bins,fs,'q',q,'atomHopFactor',atomHopFactor,'thresh',thresh,'win',winFlag,'perfRast',1);
end

%% cqt full (fast)-------------------------------------------
maxBlock = cqtKernel.fftLEN * 2^(octaveNr-1); %largest FFT Block (virtual)
suffixZeros = maxBlock; 
prefixZeros = maxBlock;
x = [zeros(prefixZeros,1); x; zeros(suffixZeros,1)]; %zeropadding
OVRLP = cqtKernel.fftLEN - cqtKernel.fftHOP;
atomNr = cqtKernel.atomNr;
bins = cqtKernel.bins;
FFTLen = cqtKernel.fftLEN;
K = cqtKernel.fKernel;
ahop = cqtKernel.atomHOP;
emptyHops = cqtKernel.firstcenter/cqtKernel.atomHOP;
    
for i=1:octaveNr  
    inc = ahop/2^(i-1);
    binVec = bins*(octaveNr-i)+1:bins*(octaveNr-i+1);
    drop = emptyHops*2^(octaveNr-i)-emptyHops; %first coefficients of all octaves have to be in synchrony
    xx = buffer(x,cqtKernel.fftLEN, OVRLP,'nodelay'); %generating FFT blocks
    XX = fft(xx); %applying fft to each column (each FFT frame)
    for n=1:2^(i-1) 
        shift = (n-1)*inc;
        phShiftVec = exp(-1i*2*pi.*(0:(FFTLen-1))'*shift./FFTLen);
        phShiftMat = repmat(phShiftVec,1,size(K,2));
        thisK = K .* phShiftMat;
        thisK = thisK';
        thisXcq = thisK*XX; %calculating cqt coefficients for all FFT frames for this octave
        if atomNr > 1
            Xoct = zeros(bins,atomNr*size(thisXcq,2)-drop);
            for u=1:bins %reshape to continous windows for each bin (for the case of several wins per frame)
               octX_bin = thisXcq((u-1)*atomNr+1:u*atomNr,:);
               Xcont = reshape(octX_bin,1,size(octX_bin,1)*size(octX_bin,2));
               Xoct(u,:) = Xcont(1+drop:end);
            end
            thisXcq = Xoct;
        else
            thisXcq = thisXcq(:,1+drop:end);
        end
        tVec = 1:2^(i-1):size(thisXcq,2)*2^(i-1);
        tVec = tVec + (n-1);
        spCQT(binVec,tVec) = thisXcq;
    end

    if i~=octaveNr
        x = filtfilt(B,A,x); %anti aliasing filter
        x = x(1:2:end); %drop samplerate by 2
    end  
end

%% return
intParam = struct('sufZeros',suffixZeros,'preZeros',prefixZeros,'xlen_init',xlen_init,'fftLEN',cqtKernel.fftLEN,'fftHOP',cqtKernel.fftHOP,...
    'q',q,'filtCoeffA',A,'filtCoeffB',B,'firstcenter',cqtKernel.firstcenter,'atomHOP',cqtKernel.atomHOP,...
    'atomNr',cqtKernel.atomNr,'Nk_max',cqtKernel.Nk_max,'Q',cqtKernel.Q,'rast',1);

Xcqt = struct('spCQT',spCQT,'fKernel',cqtKernel.fKernel,'fmax',fmax,'fmin',fmin,'octaveNr',octaveNr,...
    'bins',cqtKernel.bins,'intParams',intParam);