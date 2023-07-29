#!/bin/bash

source frame_number.sh

frame_no=$frame_no
pdb=".pdb"
file_name="sic1_frame_$frame_no$pdb"

cp "cl_2_frame_$frame_no$pdb" $file_name

shift=10000000

echo '#!/bin/bash' > resid_change.sh

echo source frame_number.sh >> resid_change.sh

echo frame_no=$frame_no >> resid_change.sh 
echo pdb="'.pdb'" >> resid_change.sh
echo file_name="'cl_1_sic1_frame_$frame_no$pdb'" >> resid_change.sh



for j in {0..7}
do 
 
res=$(( 10000000 + $j*10000000 + $shift ))

#echo $j $res 
#echo 'A   '"$j"
init_tab='A   '
final_tab='   '
str_res_init="$init_tab""$j""$final_tab"
str_res_final="$init_tab""$res""$final_tab"

echo "sed -i 's/$str_res_init/$str_res_final/g' $file_name" >> resid_change.sh

done


for j in {8..9}
do 
 
res=$(( 10000000 + $j*10000000 + $shift ))

init_tab='A   '
init_tab_1='A  '
final_tab='   '
str_res_init="$init_tab""$j""$final_tab"
str_res_final="$init_tab_1""$res""$final_tab"

echo "sed -i 's/$str_res_init/$str_res_final/g' $file_name" >> resid_change.sh

done


for j in {10..90}
do 
 
res=$(( 10000000 + $j*10000000 + $shift ))

init_tab='A  '
final_tab='   '
str_res_init="$init_tab""$j""$final_tab"
str_res_final="$init_tab""$res""$final_tab"

echo "sed -i 's/$str_res_init/$str_res_final/g' $file_name" >> resid_change.sh

done

echo "sed -i 's/A9999   /A   100   /g' $file_name" >> resid_change.sh
echo "sed -i 's/A   100/A   1/g' $file_name" >> resid_change.sh
echo "sed -i 's/0000000    /    /g' $file_name" >> resid_change.sh
echo "sed -i 's/TER    1375      THR A  90/TER    1375      THR A  92/g' $file_name" >> resid_change.sh

chmod +x resid_change.sh
./resid_change.sh

