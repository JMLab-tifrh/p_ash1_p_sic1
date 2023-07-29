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

#-----Building TOPOLOGY-----#
gmx pdb2gmx -f $filename_start$frame_no$pdb -o $filename_start$frame_no$gro -p $topol$frame_no$top -ignh < pdb2gmx_choices.txt >& log_pdb_pdb2gmx

#-----Adding the box-----#
#gmx editconf -f $filename_start$frame_no$gro -o $filename_start$frame_no$box$gro -c -d $dist_from_box_edge -bt $box_type >& log_box_editconf
gmx editconf -f $filename_start$frame_no$gro -o $filename_start$frame_no$box$gro -c -box $box_length $box_length $box_length -bt $box_type >& log_box_editconf

#-----Solvation-----#
gmx solvate -cp $filename_start$frame_no$box$gro -cs $water_model -o $filename_start$frame_no$box$solvated$gro -p $topol$frame_no$top >& log_box_solvate_solvate

#-----Neutralizing and adding NACL as per concentration in Molar units-----#
gmx grompp -f "ions.mdp" -c $filename_start$frame_no$box$solvated$gro -p $topol$frame_no$top -o "ions.tpr" -maxwarn 1 >& log_box_solvate_ions_grompp

gmx genion -s "ions.tpr" -o $filename_start$frame_no$box$solvated$ions$gro -p $topol$frame_no$top -pname NA -nname CL -neutral -conc $conc < genion_choices.txt >& log_box_solvate_ions_genion

#-----Energy minimization-----#
gmx grompp -f "minim.mdp" -c $filename_start$frame_no$box$solvated$ions$gro -p $topol$frame_no$top -o $em$tpr >& log_energy_minimize_grompp

gmx mdrun -v -deffnm $em >& log_energy_minimize_mdrun

#-----NVT equilibration-----#
gmx grompp -f "nvt.mdp" -c $em$gro -r $em$gro -p $topol$frame_no$top -o $nvt$tpr >& log_nvt_grompp

gmx mdrun -v -deffnm $nvt >& log_nvt_mdrun

#-----NPT equilibration-----#
gmx grompp -f "npt.mdp" -c $nvt$gro -r $nvt$gro -t $nvt$cpt -p $topol$frame_no$top -o $npt$tpr >& log_npt_grompp

gmx mdrun -v -deffnm $npt >& log_npt_mdrun

#-----Production MD-----#
gmx grompp -f "md.mdp" -c $npt$gro -t $npt$cpt -p $topol$frame_no$top -o $pmd$tpr >& log_pmd_grompp

gmx mdrun -v -deffnm $pmd >& log_pmd_mdrun


