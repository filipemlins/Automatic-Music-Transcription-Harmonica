

function [Y, gaps] = fill_gaps(X)


Y = X;
s = size(X);

%% number of templates per pitch
KT = 1;
if length(s)<2 | length(s)>3
    error('Template must have 2 or 3 dimensions');
elseif length(s)==3
    KT = s(3);
end

for k=1:KT;
    
    tpl = X(:,:,k);
    ss = sum(tpl);
    range = find(diff(ss)~=0);
    if length(range)<2
        range = [1 s(2)];
    end
    range = [range(1), range(end)];
    gaps = find(ss==0);
    gaps = gaps(gaps>range(1)); %% trimout 
    gaps = gaps(gaps<range(2));
   
    for i=1:length(gaps);
        prev = tpl(:,gaps(i)-1);
        tpl(:,gaps(i)) = [[0 0 0 0 0]'; prev(1:end-5)];
    end
    
    tpl = fixTemplates(tpl');     
    Y(:,:,k) = tpl';    
    
end
    




