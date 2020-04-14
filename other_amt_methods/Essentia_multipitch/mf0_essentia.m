 function [h, notes_k, notes_m] = mf0_essentia(input_wav, outfile)
 h = 0;
[path_wav,wav_name,ext] = fileparts(input_wav);

%path_essentia = '/Users/schramm/GoogleDrive/CARREIRA/UFRGS/RESEARCH/PROJECTS_UFRGS/AMT/essentia-2.1_beta3/build/src/examples/';
%path_essentia = '/Users/schramm/Desktop/automaticharmonicatranscription/tools/essentia-extractors-v2.1_beta2/';
path_essentia = 'tools/essentia-2.1_beta3/build/src/examples/';
curdir = pwd; 
cd(path_essentia);

mf0_tool = './essentia_standard_pitchdemo';
%mf0_tool = './standard_pitchdemo';
%[mf0_tool, ' ', input_wav, ' ', [outfile, '_klapuri.txt'] , '  3' ],
[status_k,result_k] = system([mf0_tool, ' ', input_wav, ' ', [outfile, '_klapuri.txt'] , '  3' ]); % klapuri==3
[status_m,result_m] = system([mf0_tool, ' ', input_wav, ' ', [outfile, '_melodia.txt'] , '  4' ]); % melodia==4

cd(curdir);
notes_k = dlmread([outfile, '_klapuri.txt']);
notes_m = dlmread([outfile, '_melodia.txt']);

h = status_k | status_m;



