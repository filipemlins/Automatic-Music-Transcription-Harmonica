function [Pre,Rec,F] = computeFrameBasedMetrics(pianoRoll,pianoRollGT_expand, pianoRollGT)


if nargin<3
    pianoRollGT = pianoRollGT_expand;
end

pianoRoll = double(pianoRoll);
pianoRollGT_expand = double(pianoRollGT_expand);

% Compute frame-based metrics
Nref = sum(sum(pianoRollGT));
Ntot = sum(sum(pianoRoll));
Ntp = sum(sum(pianoRoll+pianoRollGT_expand > 1));

Rec = Ntp/(Nref+eps);
Pre = Ntp/(Ntot+eps);
F = 2*((Pre*Rec)/(Pre+Rec+eps));