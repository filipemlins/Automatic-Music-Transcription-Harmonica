function [Pre, Rec, F, pianoRoll] = fmeasure_dynamic(plca_bin_s, gts, lim)
%

pianoRoll = zeros(size(gts));
for i=1:size(gts,2)
    pianoRoll(:,i) = double(plca_bin_s(:,i)) > lim(i);
end
    

pianoRollGT = double(gts);

% Compute frame-based metrics
Nref = sum(sum(pianoRollGT));
Ntot = sum(sum(pianoRoll));
Ntp = sum(sum(pianoRoll+pianoRollGT > 1));

Rec = Ntp/(Nref+eps);
Pre = Ntp/(Ntot+eps);
F = 2*((Pre*Rec)/(Pre+Rec+eps));
