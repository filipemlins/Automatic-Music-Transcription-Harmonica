M = dir('TrainFolderMagenta\*.*');

for i=3:(size(M,1))
filename=fullfile(M(i,1).folder,M(i,1).name);

%salvar arquivo no mesmo local com nome adicionado midi.txt
[filepath,name,ext] = fileparts(filename)

MAPS_MUS = 'MAPS_MUS-';
name_new = strcat(MAPS_MUS,name)
END_MAPS='_StbgTGd2';
name_new = strcat(name_new,END_MAPS)


filepath_new = strcat(filepath,'Renamed');

copyfile(fullfile(filepath,[name ext]),fullfile(filepath_new,[name_new ext]));

end