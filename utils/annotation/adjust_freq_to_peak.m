

function new_freq_pos = adjust_freq_to_peak(spec, freq_pos, w, bins_div)

% spec = input spectrum (cqt frequencies)
% freq_pos = pitch annotation (approximation)
% w = neighbourhood for searching
% figure
% plot(spec, 'r'); hold on;
% plot(freq_pos, spec(freq_pos), 'ro');
% plot(freq_pos+12*bins_div, spec(freq_pos+12*bins_div), 'ro');


%% smoothing the input spectrum
[bb,aa]=butter(3,.35);
spec = filtfilt(bb,aa,spec);

% plot(spec, 'b');

if freq_pos+w>size(spec,1)
    new_freq_pos=freq_pos;
    return 
end
aln = zeros(size(spec));
p = spec(freq_pos-w:freq_pos+w);
[m,idx] = max(p);

idx2=0;
if freq_pos+12*bins_div+w<=size(spec,1)
aln2 = zeros(size(spec));
p = spec((freq_pos+12*bins_div-w):(freq_pos+12*bins_div+w));
[m,idx2] = max(p);
end

%idx = round((idx+idx2)/2);

new_freq_pos = freq_pos + idx - w;

% plot(new_freq_pos+idx-w, spec(freq_pos+idx-w), 'bo');
% plot(new_freq_pos+12*bins_div+idx2-w, spec(freq_pos+12*bins_div+idx2-w), 'bo');
%figure;


