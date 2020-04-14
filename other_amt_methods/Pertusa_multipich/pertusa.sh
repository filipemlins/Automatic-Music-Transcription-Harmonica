#!/bin/bash
#
#
if [ $# -eq 1 ]
then
    sonic-annotator -t mf0ua0.n3  -w csv --csv-basedir mf0ua_output/ --csv-force  $1
else
    sonic-annotator -t mf0ua0.n3  -w csv --csv-basedir $2 --csv-force  $1
fi

