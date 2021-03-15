#/bin/bash

SA_folder="xx_x"
which_cases="x"
heart="xx"

call_meshtool="/home/common/meshtool_new/meshtool/meshtool interpolate node2elem "

directory="/media/crg17/Seagate\\ Backup\\ Plus\\ Drive/CT_cases/"$which_cases"_case"$heart
for vein in AN AL LA IL IN;
do for elec_num in `seq 1 8`
do cmd=$call_meshtool" -omsh="$directory"/meshing/1000um/BiV/BiV -idat="$directory"/simulations/multipole/init_files/"$SA_folder"_"$vein"_"$elec_num".dat -odat="$directory"/simulations/multipole/init_files/"$SA_folder"_"$vein"_"$elec_num"_elemwise.dat"
echo $cmd
eval $cmd
done
done
