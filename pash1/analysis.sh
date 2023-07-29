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

#-----Potential after energy minimization-----#
gmx energy -f $em$edr -o $potential$frame_no$xvg < potential_choices.txt >& log_energy_minimize_potential_plot

#-----Temperature after NVT equilibration-----#
gmx energy -f $nvt$edr -o $temperature$frame_no$xvg < temperature_choices.txt >& log_nvt_temperature_plot

#-----Pressure and density after NPT equilibration-----#
gmx energy -f $npt$edr -o $pressure$frame_no$xvg < pressure_choices.txt >& log_npt_pressure_plot

gmx energy -f $npt$edr -o $density$frame_no$xvg < density_choices.txt >& log_npt_density_plot

#-----PBC correction and radius of gyration after production MD-----#
gmx trjconv -pbc whole -f $pmd$xtc -o $pmd$whole$xtc -s $pmd$tpr < whole_choices.txt >& log_pmd_whole

gmx trjconv -pbc mol -center -f $pmd$whole$xtc -o $pmd$whole$nopbc$xtc -s $pmd$tpr < pbc_choices.txt >& log_pmd_whole_nopbc

gmx gyrate -s $pmd$tpr -f $pmd$whole$nopbc$xtc -o $gyrate$frame_no$xvg < gyrate_choices.txt >& log_pmd_gyration_plot



