#!/bin/bash

source frame_number.sh

#-----Input parameters-----#
gpu_id="0"
ntmpi="1"
ntomp="16"

frame_no=$frame_no
force_field="1"
water_type="1"
dist_from_box_edge="2.0"
box_length="20"
conc="0.15"
box_type="cubic"
replaced_mol="SOL"
water_model="tip4p2005.gro"

#-----naming_of_files-----#
filename_start="edited_pash1_"
pdb=".pdb"
gro=".gro"
tpr=".tpr"
top=".top"
edr=".edr"
xvg=".xvg"
cpt=".cpt"
xtc=".xtc"

topol="topol_"
box="_box"
solvated="_solv"
ions="_ions"

em="em_$frame_no"
nvt="nvt_$frame_no"
npt="npt_$frame_no"
pmd="md_$frame_no"

potential="potential_"
temperature="temperature_"
pressure="pressure_"
density="density_"
gyrate="gyrate_"
whole="_whole"
nopbc="_nopbc"

gmx mdrun -ntmpi $ntmpi -ntomp $ntomp -gpu_id $gpu_id -v -deffnm $pmd -s $pmd$tpr -cpi $pmd$cpt  >& log_pmd_mdrun

