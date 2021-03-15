#/bin/bash

SA_folder="xx_x"
which_cases="x"
heart="xx"
midseptum="false"

call_meshtool="/home/common/meshtool_new/meshtool/meshtool interpolate node2elem "

directory="/media/crg17/Seagate\\ Backup\\ Plus\\ Drive/CT_cases/"$which_cases"_case"$heart

if [[ "$midseptum" == "false" ]]
then
cmd=$call_meshtool" -omsh="$directory"/meshing/1000um/BiV/BiV -idat="$directory"/simulations/multipole/init_files/"$SA_folder"_BiV.RV_endo.apex.dat -odat="$directory"/simulations/multipole/init_files/"$SA_folder"_BiV.RV_endo.apex_elemwise.dat"
elif [[ "$midseptum" == "true" ]]
then
cmd=$call_meshtool" -omsh="$directory"/meshing/1000um/BiV/BiV -idat="$directory"/simulations/multipole/init_files/"$SA_folder"_BiV.midseptum.dat -odat="$directory"/simulations/multipole/init_files/"$SA_folder"_BiV.midseptum_elemwise.dat"
fi
echo $cmd
eval $cmd
