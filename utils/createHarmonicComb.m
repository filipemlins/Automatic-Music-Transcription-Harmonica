function [harmonicComb] = createHarmonicComb()

% With 60 bins/octave log-frequency resolution...

% Initialize
harmonicComb = zeros(600,96);
combTemplate = round(60*log2(1:30));

% For each note
for i=60:96
    
    f0_index = i*5 - 21*5 ;
    
    comb1 = f0_index + (combTemplate-2);
    comb2 = f0_index + (combTemplate-1);
    comb3 = f0_index + (combTemplate);
    comb4 = f0_index + (combTemplate+1);
    comb5 = f0_index + (combTemplate+2);
    comb6 = f0_index + (combTemplate+3);
    comb7 = f0_index + (combTemplate-3);
    comb8 = f0_index + (combTemplate+4);
    comb9 = f0_index + (combTemplate-4);    
    comb10 = f0_index + (combTemplate+5);
    comb11 = f0_index + (combTemplate-5);     
    
    tempComb = [comb1 comb2 comb3 comb4 comb5 comb6 comb7 comb8 comb9 comb10 comb11];
    
    activeIndices = find(tempComb < 546);
    
    harmonicComb(tempComb(activeIndices),i) = 1;
    
end

harmonicComb = harmonicComb(1:545,:);