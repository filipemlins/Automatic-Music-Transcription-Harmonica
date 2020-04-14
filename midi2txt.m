

addpath('MagentaExperiment/MagentaResults/');
M = dir('**/*.mid');
delimiterIn=',';

for i=1:size(M,1)
    filename=fullfile(M(i,1).folder,M(i,1).name)
    midi = readmidi(filename);

    midinotes = midiInfo(midi,0);
    [m,n] = size(midinotes);
    array = zeros(m,3);

    array(:,3) = midinotes(:,3);
    array(:,1) = midinotes(:,5);
    array(:,2) = midinotes(:,6);
    
 
    for w=1:m 
        array(w,3) = midi2freq(array(w,3));
    end
    dlmwrite(fullfile(filepath,[name,'.out']),A,delimiterIn);
end


