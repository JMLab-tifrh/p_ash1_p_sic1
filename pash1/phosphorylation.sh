#!/bin/bash

source frame_number.sh

frame_no=$frame_no
pdb=".pdb"
file_name="edited_pash1_$frame_no$pdb"

cp "pash1_$frame_no$pdb" $file_name

sed -i 's/SER P   7/S1P P   7/g' $file_name
sed -i 's/SER P   9/S1P P   9/g' $file_name
sed -i 's/THR P  12/T1P P  12/g' $file_name
sed -i 's/SER P  25/S1P P  25/g' $file_name
sed -i 's/THR P  33/T1P P  33/g' $file_name
sed -i 's/SER P  35/S1P P  35/g' $file_name
sed -i 's/SER P  38/S1P P  38/g' $file_name
sed -i 's/SER P  48/S1P P  48/g' $file_name
sed -i 's/SER P  52/S1P P  52/g' $file_name
sed -i 's/SER P  73/S1P P  73/g' $file_name
sed -i 's/HSD/HID/g' $file_name

