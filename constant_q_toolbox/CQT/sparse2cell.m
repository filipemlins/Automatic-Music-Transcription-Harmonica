function cellCQT = sparse2cell(spCQT,bins,octaveNr,atomNr,firstcenter,atomHOP)
% Maps the sparse matrix respresentation of the CQT coefficients back to
% the cell representation for inverse transform
%
%Christian Schörkhuber, Anssi Klapuri 2010-06

emptyHops = firstcenter/atomHOP; %void atom hopsizes in the beginning of the temporal kernel
cellCQT = cell(octaveNr,1);

for i=1:octaveNr
    dropped = emptyHops*2^(octaveNr-i)-emptyHops;
    X = full(spCQT(bins*octaveNr-i*bins+1:bins*octaveNr-(i-1)*bins,1:2^(i-1):end));
    X = [zeros(bins,dropped) X]; 
    X = [X zeros(bins,ceil(size(X,2)/atomNr)*atomNr-size(X,2))];
    if atomNr > 1 %reshape
        Xcell = zeros(bins*atomNr,ceil(size(X,2)/atomNr));    
        for u=1:bins;  
            Xbin = reshape(X(u,:),atomNr,length(X(u,:))/atomNr);
            Xcell((u-1)*atomNr+1:u*atomNr,:) = Xbin;
        end      
        cellCQT{i} = Xcell;
    else
        cellCQT{i} = full(X);
    end
end


