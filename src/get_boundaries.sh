#!/bin/bash

# run as bash get_boundaries.sh input output
strand=$3

grep -P ".*\tCDS\t[0-9]*\t[0-9]*\t\.\t\\$strand\t.*" $1 | awk '{print $1, $4, $5}' > $2
