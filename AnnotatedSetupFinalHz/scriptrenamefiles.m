M = dir('**/*.txt');
delimiterIn=',';
headerlinesIn=0;

for i=1:size(M,1)
filename=fullfile(M(i,1).folder,M(i,1).name);

A = importdata(filename, delimiterIn, headerlinesIn);

%salvar arquivo no mesmo local com nome adicionado midi.txt
[filepath,name,ext] = fileparts(filename);
dlmwrite(fullfile(filepath,[name,'Hz.txt']),A,delimiterIn);

end