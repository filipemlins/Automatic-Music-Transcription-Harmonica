function sparKernel= sparseKernel(minFreq, maxFreq, bins, fs, thresh)

if nargin<5 thresh= 0.0054; end % for Hamming window

Q= 1/(2^(1/bins)-1); %(3)
K= ceil( bins * log2(maxFreq/minFreq) ); %(2)
fftLen= 2^nextpow2( ceil(Q*fs/minFreq) );
tempKernel= zeros(fftLen, 1);
sparKernel= [];
for k= K:-1:1;
    len= ceil( Q * fs / (minFreq*2^((k-1)/bins)) ); %(4)
    tempKernel(1:len)= hamming(len)/len .* exp(2*pi*i*Q*(0:len-1)'/len); %(7)
    specKernel= fft(tempKernel); %(8)
    specKernel(find(abs(specKernel)<=thresh))= 0; %% important tunning option
    sparKernel= sparse([specKernel sparKernel]);
end

%% normalizing
sparKernel= conj(sparKernel) / fftLen; %(10).1
