function hz = midi2freq(midi)

hz =  440*2.^((midi-69)/12);