#!/bin/bash

for i in $(find . -type f  -iname  "*.wav")
do
    /home/filipemlins/Desktop/madmom/bin/PianoTranscriptor --mirex single $i > ${i}.out  
	
done
