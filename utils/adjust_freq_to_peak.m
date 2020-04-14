

function [new_freq_pos, shift_step] = adjust_freq_to_peak(spec, freq_pos, w)

% spec = input spectrum (cqt frequencies)
% freq_pos = pitch annotation (approximation)
% w = searching neighbourhood 

DEBUG = 0;

if DEBUG
    plot(spec, 'r'); hold on;
    plot(freq_pos, spec(freq_pos), 'ro');
end

%% smoothing the input spectrum
[bb,aa]=butter(3,.35);
spec = filtfilt(bb,aa,spec);

aln = zeros(size(spec));
p = spec(freq_pos-w:freq_pos+w);
[m,idx] = max(p);
new_freq_pos = freq_pos + idx - (w+1);
shift_step = idx - (w+1);

if DEBUG
    plot(spec, 'b');
    plot(new_freq_pos, spec(new_freq_pos), 'bo');
end


