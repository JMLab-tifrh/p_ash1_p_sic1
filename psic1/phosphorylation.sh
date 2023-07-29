#!/bin/bash

#./resid.sh

source frame_number.sh

frame_no=$frame_no
pdb=".pdb"
file_name="edited_psic1_$frame_no$pdb"

cp "psic1_$frame_no$pdb" $file_name

sed -i 's/THR P   7/TPO P   7/g' $file_name
sed -i 's/THR P  35/TPO P  35/g' $file_name
sed -i 's/THR P  47/TPO P  47/g' $file_name
sed -i 's/SER P  71/SEP P  71/g' $file_name
sed -i 's/SER P  78/SEP P  78/g' $file_name
sed -i 's/SER P  82/SEP P  82/g' $file_name


