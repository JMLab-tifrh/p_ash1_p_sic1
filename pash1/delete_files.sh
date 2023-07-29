[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"
#!/bin/bash

source frame_number.sh

frame_no=$frame_no
gro=".gro"
pdb=".pdb"
filename_start="edited_pash1_"
box="_box"
solvated="_solv"
ions="_ions"
topol="topol_"
top=".top"
md="md_"
npt="npt_"
nvt="nvt_"
em="em_"

rm *.tpr
rm *.trr
rm *.edr
rm *.cpt
rm *.xtc
rm *.log
rm *#
rm *mdout.mdp
rm posre.itp
rm $filename_start$frame_no$pdb
rm $filename_start$frame_no$gro
rm $filename_start$frame_no$box$gro
rm $filename_start$frame_no$box$solvated$gro
rm $filename_start$frame_no$box$solvated$ions$gro
rm $md$frame_no$gro
rm $nvt$frame_no$gro
rm $npt$frame_no$gro
rm $em$frame_no$gro
rm *.log
rm $topol$frame_no$top
rm log_*
rm *.xvg
rm *step*
rm *.o*
rm *.e*
