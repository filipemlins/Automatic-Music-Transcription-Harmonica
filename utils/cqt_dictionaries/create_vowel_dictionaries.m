% function create_vowel_dictionaries()
%
% This function create dictionary templates based onVQT/CQT spectrogram from wav audio files. 
% input  data (wav audio files)
% output data (mat files saved on current directory).
%
%  See also: extract_cqt_from_wav_dataset, cqt_clustering, create_vowel_dictionaries
%
% Created by Rodrigo Schramm on 27/09/2016.


p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)    
    parpool(12);
end

%% soprano
[cqt_soprano_a, d_soprano_a] = extract_cqt_from_wav_dataset([], 'soprano', 'a', 'm', 'n');
save('cqt_soprano_a','cqt_soprano_a','d_soprano_a');
clear('cqt_soprano_a','cqt_soprano_a');

[cqt_soprano_e, d_soprano_e] = extract_cqt_from_wav_dataset([], 'soprano', 'e', 'm', 'n');
save('cqt_soprano_e','cqt_soprano_e','d_soprano_e');
clear('cqt_soprano_e','cqt_soprano_e');

[cqt_soprano_i, d_soprano_i] = extract_cqt_from_wav_dataset([], 'soprano', 'i', 'm', 'n');
save('cqt_soprano_i','cqt_soprano_i','d_soprano_i');
clear('cqt_soprano_i','cqt_soprano_i');

[cqt_soprano_o, d_soprano_o] = extract_cqt_from_wav_dataset([], 'soprano', 'o', 'm', 'n');
save('cqt_soprano_o','cqt_soprano_o','d_soprano_o');
clear('cqt_soprano_o','cqt_soprano_o');

[cqt_soprano_u, d_soprano_u] = extract_cqt_from_wav_dataset([], 'soprano', 'u', 'm', 'n');
save('cqt_soprano_u','cqt_soprano_u','d_soprano_u');
clear('cqt_soprano_u','cqt_soprano_u');

%% alto
[cqt_alto_a, d_alto_a] = extract_cqt_from_wav_dataset([], 'alto', 'a', 'm', 'n');
save('cqt_alto_a','cqt_alto_a','d_alto_a');
clear('cqt_alto_a','cqt_alto_a');

[cqt_alto_e, d_alto_e] = extract_cqt_from_wav_dataset([], 'alto', 'e', 'm', 'n');
save('cqt_alto_e','cqt_alto_e','d_alto_e');
clear('cqt_alto_e','cqt_alto_e');

[cqt_alto_i, d_alto_i] = extract_cqt_from_wav_dataset([], 'alto', 'i', 'm', 'n');
save('cqt_alto_i','cqt_alto_i','d_alto_i');
clear('cqt_alto_i','cqt_alto_i');

[cqt_alto_o, d_alto_o] = extract_cqt_from_wav_dataset([], 'alto', 'o', 'm', 'n');
save('cqt_alto_o','cqt_alto_o','d_alto_o');
clear('cqt_alto_o','cqt_alto_o');

[cqt_alto_u, d_alto_u] = extract_cqt_from_wav_dataset([], 'alto', 'u', 'm', 'n');
save('cqt_alto_u','cqt_alto_u','d_alto_u');
clear('cqt_alto_u','cqt_alto_u');

%% tenor
[cqt_tenor_a, d_tenor_a] = extract_cqt_from_wav_dataset([], 'tenor', 'a', 'm', 'n');
save('cqt_tenor_a','cqt_tenor_a','d_tenor_a');
clear('cqt_tenor_a','cqt_tenor_a');

[cqt_tenor_e, d_tenor_e] = extract_cqt_from_wav_dataset([], 'tenor', 'e', 'm', 'n');
save('cqt_tenor_e','cqt_tenor_e','d_tenor_e');
clear('cqt_tenor_e','cqt_tenor_e');

[cqt_tenor_i, d_tenor_i] = extract_cqt_from_wav_dataset([], 'tenor', 'i', 'm', 'n');
save('cqt_tenor_i','cqt_tenor_i','d_tenor_i');
clear('cqt_tenor_i','cqt_tenor_i');

[cqt_tenor_o, d_tenor_o] = extract_cqt_from_wav_dataset([], 'tenor', 'o', 'm', 'n');
save('cqt_tenor_o','cqt_tenor_o','d_tenor_o');
clear('cqt_tenor_o','cqt_tenor_o');

[cqt_tenor_u, d_tenor_u] = extract_cqt_from_wav_dataset([], 'tenor', 'u', 'm', 'n');
save('cqt_tenor_u','cqt_tenor_u','d_tenor_u');
clear('cqt_tenor_u','cqt_tenor_u');

%% baritone
[cqt_baritone_a, d_baritone_a] = extract_cqt_from_wav_dataset([], 'baritone', 'a', 'm', 'n');
save('cqt_baritone_a','cqt_baritone_a','d_baritone_a');
clear('cqt_baritone_a','cqt_baritone_a');

[cqt_baritone_e, d_baritone_e] = extract_cqt_from_wav_dataset([], 'baritone', 'e', 'm', 'n');
save('cqt_baritone_e','cqt_baritone_e','d_baritone_e');
clear('cqt_baritone_e','cqt_baritone_e');

[cqt_baritone_i, d_baritone_i] = extract_cqt_from_wav_dataset([], 'baritone', 'i', 'm', 'n');
save('cqt_baritone_i','cqt_baritone_i','d_baritone_i');
clear('cqt_baritone_i','cqt_baritone_i');

[cqt_baritone_o, d_baritone_o] = extract_cqt_from_wav_dataset([], 'baritone', 'o', 'm', 'n');
save('cqt_baritone_o','cqt_baritone_o','d_baritone_o');
clear('cqt_baritone_o','cqt_baritone_o');

[cqt_baritone_u, d_baritone_u] = extract_cqt_from_wav_dataset([], 'baritone', 'u', 'm', 'n');
save('cqt_baritone_u','cqt_baritone_u','d_baritone_u');
clear('cqt_baritone_u','cqt_baritone_u');

%% bass
[cqt_bass_a, d_bass_a] = extract_cqt_from_wav_dataset([], 'bass', 'a', 'm', 'n');
save('cqt_bass_a','cqt_bass_a','d_bass_a');
clear('cqt_bass_a','cqt_bass_a');

[cqt_bass_e, d_bass_e] = extract_cqt_from_wav_dataset([], 'bass', 'e', 'm', 'n');
save('cqt_bass_e','cqt_bass_e','d_bass_e');
clear('cqt_bass_e','cqt_bass_e');

[cqt_bass_i, d_bass_i] = extract_cqt_from_wav_dataset([], 'bass', 'i', 'm', 'n');
save('cqt_bass_i','cqt_bass_i','d_bass_i');
clear('cqt_bass_i','cqt_bass_i');

[cqt_bass_o, d_bass_o] = extract_cqt_from_wav_dataset([], 'bass', 'o', 'm', 'n');
save('cqt_bass_o','cqt_bass_o','d_bass_o');
clear('cqt_bass_o','cqt_bass_o');

[cqt_bass_u, d_bass_u] = extract_cqt_from_wav_dataset([], 'bass', 'u', 'm', 'n');
save('cqt_bass_u','cqt_bass_u','d_bass_u');
clear('cqt_bass_u','cqt_bass_u');


%%===============================================================
%%===============================================================


%% -- generate dictionaries
% SOPRANO
load cqt_soprano_a; 
[dic_soprano_a, err_soprano_a] = cqt_clustering(cqt_soprano_a, d_soprano_a); 
save('dic_soprano_a', 'dic_soprano_a', 'err_soprano_a');
clear('dic_soprano_a', 'd_soprano_a', 'err_soprano_a');

load cqt_soprano_e; 
[dic_soprano_e, err_soprano_e] = cqt_clustering(cqt_soprano_e, d_soprano_e); 
save('dic_soprano_e', 'dic_soprano_e', 'err_soprano_e');
clear('dic_soprano_e', 'd_soprano_e', 'err_soprano_e');

load cqt_soprano_i; 
[dic_soprano_i, err_soprano_i] = cqt_clustering(cqt_soprano_i, d_soprano_i); 
save('dic_soprano_i', 'dic_soprano_i', 'err_soprano_i');
clear('dic_soprano_i', 'd_soprano_i', 'err_soprano_i');

load cqt_soprano_o; 
[dic_soprano_o, err_soprano_o] = cqt_clustering(cqt_soprano_o, d_soprano_o); 
save('dic_soprano_o', 'dic_soprano_o', 'err_soprano_o');
clear('dic_soprano_o', 'd_soprano_o', 'err_soprano_o');

load cqt_soprano_u; 
[dic_soprano_u, err_soprano_u] = cqt_clustering(cqt_soprano_u, d_soprano_u); 
save('dic_soprano_u', 'dic_soprano_u', 'err_soprano_u');
clear('dic_soprano_u', 'd_soprano_u', 'err_soprano_u');

%% ALTO
load cqt_alto_a; 
[dic_alto_a, err_alto_a] = cqt_clustering(cqt_alto_a, d_alto_a); 
save('dic_alto_a', 'dic_alto_a', 'err_alto_a');
clear('dic_alto_a', 'd_alto_a', 'err_alto_a');

load cqt_alto_e; 
[dic_alto_e, err_alto_e] = cqt_clustering(cqt_alto_e, d_alto_e); 
save('dic_alto_e', 'dic_alto_e', 'err_alto_e');
clear('dic_alto_e', 'd_alto_e', 'err_alto_e');

load cqt_alto_i; 
[dic_alto_i, err_alto_i] = cqt_clustering(cqt_alto_i, d_alto_i); 
save('dic_alto_i', 'dic_alto_i', 'err_alto_i');
clear('dic_alto_i', 'd_alto_i', 'err_alto_i');

load cqt_alto_o; 
[dic_alto_o, err_alto_o] = cqt_clustering(cqt_alto_o, d_alto_o); 
save('dic_alto_o', 'dic_alto_o', 'err_alto_o');
clear('dic_alto_o', 'd_alto_o', 'err_alto_o');

load cqt_alto_u; 
[dic_alto_u, err_alto_u] = cqt_clustering(cqt_alto_u, d_alto_u); 
save('dic_alto_u', 'dic_alto_u', 'err_alto_u');
clear('dic_alto_u', 'd_alto_u', 'err_alto_u');

%% TENOR
load cqt_tenor_a; 
[dic_tenor_a, err_tenor_a] = cqt_clustering(cqt_tenor_a, d_tenor_a); 
save('dic_tenor_a', 'dic_tenor_a', 'err_tenor_a');
clear('dic_tenor_a', 'd_tenor_a', 'err_tenor_a');

load cqt_tenor_e; 
[dic_tenor_e, err_tenor_e] = cqt_clustering(cqt_tenor_e, d_tenor_e); 
save('dic_tenor_e', 'dic_tenor_e', 'err_tenor_e');
clear('dic_tenor_e', 'd_tenor_e', 'err_tenor_e');

load cqt_tenor_i; 
[dic_tenor_i, err_tenor_i] = cqt_clustering(cqt_tenor_i, d_tenor_i); 
save('dic_tenor_i', 'dic_tenor_i', 'err_tenor_i');
clear('dic_tenor_i', 'd_tenor_i', 'err_tenor_i');

load cqt_tenor_o; 
[dic_tenor_o, err_tenor_o] = cqt_clustering(cqt_tenor_o, d_tenor_o); 
save('dic_tenor_o', 'dic_tenor_o', 'err_tenor_o');
clear('dic_tenor_o', 'd_tenor_o', 'err_tenor_o');

load cqt_tenor_u; 
[dic_tenor_u, err_tenor_u] = cqt_clustering(cqt_tenor_u, d_tenor_u); 
save('dic_tenor_u', 'dic_tenor_u', 'err_tenor_u');
clear('dic_tenor_u', 'd_tenor_u', 'err_tenor_u');

%% BARITONE
load cqt_baritone_a; 
[dic_baritone_a, err_baritone_a] = cqt_clustering(cqt_baritone_a, d_baritone_a); 
save('dic_baritone_a', 'dic_baritone_a', 'err_baritone_a');
clear('dic_baritone_a', 'd_baritone_a', 'err_baritone_a');

load cqt_baritone_e; 
[dic_baritone_e, err_baritone_e] = cqt_clustering(cqt_baritone_e, d_baritone_e); 
save('dic_baritone_e', 'dic_baritone_e', 'err_baritone_e');
clear('dic_baritone_e', 'd_baritone_e', 'err_baritone_e');

load cqt_baritone_i; 
[dic_baritone_i, err_baritone_i] = cqt_clustering(cqt_baritone_i, d_baritone_i); 
save('dic_baritone_i', 'dic_baritone_i', 'err_baritone_i');
clear('dic_baritone_i', 'd_baritone_i', 'err_baritone_i');

load cqt_baritone_o; 
[dic_baritone_o, err_baritone_o] = cqt_clustering(cqt_baritone_o, d_baritone_o); 
save('dic_baritone_o', 'dic_baritone_o', 'err_baritone_o');
clear('dic_baritone_o', 'd_baritone_o', 'err_baritone_o');

load cqt_baritone_u; 
[dic_baritone_u, err_baritone_u] = cqt_clustering(cqt_baritone_u, d_baritone_u); 
save('dic_baritone_u', 'dic_baritone_u', 'err_baritone_u');
clear('dic_baritone_u', 'd_baritone_u', 'err_baritone_u');

%% BASS
load cqt_bass_a; 
[dic_bass_a, err_bass_a] = cqt_clustering(cqt_bass_a, d_bass_a); 
save('dic_bass_a', 'dic_bass_a', 'err_bass_a');
clear('dic_bass_a', 'd_bass_a', 'err_bass_a');

load cqt_bass_e; 
[dic_bass_e, err_bass_e] = cqt_clustering(cqt_bass_e, d_bass_e); 
save('dic_bass_e', 'dic_bass_e', 'err_bass_e');
clear('dic_bass_e', 'd_bass_e', 'err_bass_e');

load cqt_bass_i; 
[dic_bass_i, err_bass_i] = cqt_clustering(cqt_bass_i, d_bass_i); 
save('dic_bass_i', 'dic_bass_i', 'err_bass_i');
clear('dic_bass_i', 'd_bass_i', 'err_bass_i');

load cqt_bass_o; 
[dic_bass_o, err_bass_o] = cqt_clustering(cqt_bass_o, d_bass_o); 
save('dic_bass_o', 'dic_bass_o', 'err_bass_o');
clear('dic_bass_o', 'd_bass_o', 'err_bass_o');

load cqt_bass_u; 
[dic_bass_u, err_bass_u] = cqt_clustering(cqt_bass_u, d_bass_u); 
save('dic_bass_u', 'dic_bass_u', 'err_bass_u');
clear('dic_bass_u', 'd_bass_u', 'err_bass_u');


% remove parpool job
poolobj = gcp('nocreate');
delete(poolobj);




