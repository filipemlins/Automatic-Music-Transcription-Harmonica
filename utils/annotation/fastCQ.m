
function cq = fastCQ(x, sparKernel) % x must be a row vector
%cq= fft(x,size(sparKernel,1)) * sparKernel;
fr = fft(x,size(sparKernel,1));
cq= fr * sparKernel;
