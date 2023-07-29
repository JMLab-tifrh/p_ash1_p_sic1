#!/bin/bash

source frame_number.sh

frame_no=$frame_no

gmx mindist -pi yes -f md_"$frame_no"_whole_nopbc.xtc  -s md_"$frame_no".tpr -od min_dist_"$frame_no".xvg<<EOF
1
EOF
