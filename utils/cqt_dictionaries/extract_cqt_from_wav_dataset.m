function [cqt_mat, d] = extract_cqt_from_wav_dataset(path_wav, voice_type, vowel, dynamic, var_type)
%
% This function extracts VQT/CQT spectrogram from wav audio files. 
%  Inputs:: 
%     Please refer to function rwc_path_filter (help rwc_path_filter).
%  Outputs:: 
%     cqt_mat: a cell array which holds the VQT/CQT spectrograms.
%     d: the respective wav file path used on the spectrogram extraction.
%
%  See also: cqt_clustering, create_vowel_dictionaries
%
% Created by Rodrigo Schramm on 27/09/2016.

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)    
    parpool(12);
end



path_wav = '/import/c4dm-04/schramm/RWA_dataset/';
%path_wav = '/Users/schramm/Desktop/RWA_dataset/';
[d, ~] = rwc_path_filter(path_wav, voice_type, vowel, dynamic, var_type);

max_files = 100;
cqt_mat = {};

if length(d)>max_files
   warning('Number of files is exceeding max_files (%d > %d). parfoor loop will be limited to: %d', length(d),  max_files, max_files);
end

 parfor i=1:min(length(d), max_files)   
   wav_file = d(i).name;
   fprintf('file: %s\n', wav_file);
    
   if 0
    % Compute CQT and perform simple noise reduction
    intCQT = computeCQT(wav_file);    
    X = intCQT(:,round(1:7.1128:size(intCQT,2)))';
    noiseLevel1 = medfilt1(X',40);
    noiseLevel2 = medfilt1(min(X',noiseLevel1),40);
    X = max(X-noiseLevel2',0);
    Y = X(1:4:size(X,1),:);
    
   else % VQT has better results on low frequencies (less blur on F0)
    
    % Compute VQT and perform simple noise reduction
    [intVQT] = computeVQT(wav_file);
    X = intVQT(:,round(1:5.2883:size(intVQT,2)))';
    noiseLevel1 = medfilt1(X',40);
    noiseLevel2 = medfilt1(min(X',noiseLevel1),40);
    X = max(X-noiseLevel2',0);
    Y = X(1:2:size(X,1),:);    
   end
       
    %clear('intCQT','X','noiseLevel1','noiseLevel2');
    
    cqt_mat{i} = permute(Y, [2 1 3]);  % 40ms step    
 end


% remove parpool job
%poolobj = gcp('nocreate');
%delete(poolobj);






