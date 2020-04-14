function txt2midi(filename)

delimiterIn=',';
headerlinesIn=0;
A = importdata(filename, delimiterIn, headerlinesIn);

M = zeros(size(A,1),6); %cria matriz de zeros

M(:,1) = 1;         % all in track 1
M(:,2) = 1;         % all in channel 1
M(:,4) = 60;   %volume


for i=1:size(A,1)
M(i,3) = freq2midi(A(i,3)); % random note numbers
M(i,5) = A(i,1);
M(i,6) = A(i,2);
end

midi_new = matrix2midi(M);
[filepath,name,ext] = fileparts(filename);
writemidi(midi_new, fullfile(filepath,[name,'.mid']));

end