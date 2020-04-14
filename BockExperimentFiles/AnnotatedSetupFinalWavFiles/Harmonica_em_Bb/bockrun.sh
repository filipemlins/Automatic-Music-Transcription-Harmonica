#!/bin/bash

source activate magenta
for i in $(find . -type f  -iname  "*.wav")
do
    /home/filipemlins/Desktop/madmom/bin/PianoTranscriptor --mirex single $i > ${i}.out  
    /home/filipemlins/Desktop/madmom/bin/PianoTranscriptor --midi single $i > ${i}.mid  	
done
