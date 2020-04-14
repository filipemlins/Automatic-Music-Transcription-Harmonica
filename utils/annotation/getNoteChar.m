function c = getNoteChar(midi_pitch)

%notes = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
notes = {'C','Db','D','Eb','E','F','Gb','G','Ab','A','Bb','B'};
idx = mod(midi_pitch,12)+1;
c = notes{idx};


