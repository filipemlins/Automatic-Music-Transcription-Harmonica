
M = dir('**/*Hz.txt');
delimiterIn=',';
headerlinesIn=0;

for i=1:size(M,1)

filename=fullfile(M(i,1).folder,M(i,1).name);

A = importdata(filename, delimiterIn, headerlinesIn);
size(A,1);
B = A;
for j=1:size(A,1)
B(j,3) = round(freq2midi(A(j,3)));
end
%salvar arquivo no mesmo local com nome adicionado midi.txt

[filepath,name,ext] = fileparts(filename);
name= strrep(name,'Hz','')

file=fullfile(filepath,[name,'midi.txt']);
FID = fopen(file,'w');
if (FID == -1)
    error('Cannot open file %s', fileName)
end
fprintf(FID, "Onset Time\tOffset Time\tMidi Pitch\n");

for z=1:size(B,1)
   
    fprintf(FID,'%.2f\t%.2f\t%g\n',B(z,1), B(z,2), B(z,3));
    
end

fclose(FID);

end