#!/bin/bash

# run as bash get_boundaries.sh input output

grep -P ".*\tCDS\t[0-9]*\t[0-9]*\t\.\t\+\t.*" $1 | awk '{print $1, $4, $5}' > $2