function intCQT = schramm_cqt(y,sparKernel)
win_size=4096;
hop_size=512;
BUFFER_LEN = length(y);

c=0;
ff=[];
p=win_size+1;
while p<BUFFER_LEN;
    c=c+1;    
    ff(:,c)=fastCQ(y(p-win_size:p)',sparKernel);    
    p=p+hop_size;
end
%intCQT = abs(ff).*2^15; %% fftw do not normalise, 
% so it is need to divide by n in C++ code
intCQT = abs(ff);  

