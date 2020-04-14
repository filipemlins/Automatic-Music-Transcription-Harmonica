addpath('D:\grive\Doutorado\Ismir 2018\VersaoFinal15082018\amt_harmonica_ismir20182\utils\matlab-midi-master\src');
M = dir('**/*.txt');

for i=1:size(M,1)
    filename=fullfile(M(i,1).folder,M(i,1).name)
    txt2midi(filename);
end