function y = icqt(Xcqt)
%y = icqt(Xcqt) computes the inverse CQT of the CQT coefficients in Xcqt.spCQT
%
%The input structue Xcqt is the structure gained by cqt() and cqtPerfectRast(), respectively. 
%If the CQT coefficients in Xcqt.spCQT are not changed, the output y is the
%reconstructed (near-perfect) time-domain signal of the input signal x
%(cqt(x,...)) withing the frequency range [fmin fmax].
%
%Christian Schörkhuber, Anssi Klapuri 2010-06

cellCQT = sparse2cell(Xcqt.spCQT,Xcqt.bins,Xcqt.octaveNr,Xcqt.intParams.atomNr,Xcqt.intParams.firstcenter,Xcqt.intParams.atomHOP);

FFTLen = Xcqt.intParams.fftLEN;
octaveNr = Xcqt.octaveNr;
HOPSZ = Xcqt.intParams.fftHOP;

%% Kernel for inverse transform
Kinv = Xcqt.fKernel;

%% inverse transform
y = [];
for i = octaveNr:-1:1
    cellCQT_oct = cellCQT{i};    
    Y = Kinv * cellCQT_oct; %compute spectrum of reconstructed signal for all coefficients in this octave  
    y_oct_temp = ifft(Y);
    y_oct = 2*real(y_oct_temp); %Y contains no negative frequencies -> keep only real part*2 to 
                                %reconstruct real valued time signal 
    NBLOCKS = size(Y,2);      
    siglen = FFTLen + (NBLOCKS-1)*HOPSZ;
    y = [y;zeros(siglen-length(y),1)];
    for n = 1:NBLOCKS
        y((n-1)*HOPSZ+1:((n-1)*HOPSZ)+FFTLen) = y_oct(:,n) + y((n-1)*HOPSZ+1:((n-1)*HOPSZ)+FFTLen); %overlap-add
    end
    
    if(i~=1) %upsampling by factor two
         y = upsample(y,2); %insert one zero between each sample
         y = filtfilt(Xcqt.intParams.filtCoeffB,Xcqt.intParams.filtCoeffA,y);
         y = y * 2;
    end

end

y = y(Xcqt.intParams.preZeros+1:end); %crop introduced zeros at the beginning
y = y(1:Xcqt.intParams.xlen_init); %crop overhead zeros at the end
