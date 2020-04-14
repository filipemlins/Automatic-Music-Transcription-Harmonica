function [Pre, Rec, F] = compute_fmeasure(zz, gtzz, thr)

[Pre,Rec,F] = computeFrameBasedMetrics(zz>thr, gtzz>0);
