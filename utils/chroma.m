function c = chroma(spec)

%s = length(spec);
s = size(spec);
idx = 1:60:s(1);
c=zeros([12*5 s(2)]);



for i=1:(12*5);
    d =(i-1)+idx;
    d = d(d<s(1));
    if length(d)>0;        
        c(i,:) = sum(spec(d,:));
    end
end
    