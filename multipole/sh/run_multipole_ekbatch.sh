#/bin/bash

clear

which_cases="both"
biv_only=false
midseptum=false
all_cases="01 02 03 04 05 06 07 08 09 "`seq 10 20`
SA_folder="default_HF"
FEC_lvl="33"
specific_case=""
suffix="noPVTV"

function Usage(){
    echo "Script to run the multipole simulations using ekbatch and the init files."
    echo "Parameters:"
    echo "--specific_case: Optional argument to run only that case."
    echo "--FEC_lvl: Height of the FEC folder"
    echo "--SA_folder: Name/preffix of the simulations"
    echo "--which_cases: h/RR, HF or both for healthy/reverse remodelled, heart failure of both, respectively."
    echo "--biv_only: If true, only created for the BiV RV apex, false creates all the others. Default is false."
    echo "--midseptum: If true, it creates the init file for the midseptum electrode."
    echo "--suffix: String with the suffix of the mesh in the meshes folder."
    echo "-h/--help: Parameters usage."
}


while [ "$1" != "" ]; do
    case $1 in
        --SA_folder)            shift
                                SA_folder=$1
                                ;;
        --FEC_lvl)              shift
                                FEC_lvl=$1
                                ;;
        --specific_case)        shift
                                specific_case=$1
                                ;;
        --which_cases)          shift
                                which_cases=$1
                                ;;
        --biv_only )            shift
                                biv_only=$1
                                ;;
        --suffix )              shift
                                suffix=$1
                                ;;
        --midseptum )           shift
                                midseptum=$1
                                ;;
        -h | --help )           Usage
                                exit
                                ;;
        * )                     echo "Command not found. Use -h or --help for more info."
                                exit 1
    esac
    shift
done

if [[ $which_cases == "both" ]];
then
cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/run_multipole_ekbatch.sh --which_cases RR --biv_only "$biv_only" --midseptum "$midseptum" --SA_folder "$SA_folder" --FEC_lvl "$FEC_lvl" --suffix "$suffix
echo $cmd
eval $cmd

cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/run_multipole_ekbatch.sh --which_cases HF --biv_only "$biv_only" --midseptum "$midseptum" --SA_folder "$SA_folder" --FEC_lvl "$FEC_lvl" --suffix "$suffix
echo $cmd
eval $cmd

exit

elif [[ $which_cases  == "RR" ]]
then
which_cases="h"

elif [[ $which_cases == "HF" ]]
then
all_cases=$all_cases" 21 22 23 24"

fi

if [[ $specific_case -ne "" ]]
then
all_cases=$specific_case
fi

common_path="/media/crg17/Seagate\ Backup\ Plus\ Drive/CT_cases"

for heart in $all_cases
    do
    init_files_path=$common_path"/"$which_cases"_case"$heart"/simulations/multipole/init_files"
    mesh_path=$common_path"/"$which_cases"_case"$heart"/meshing/1000um/BiV/meshes"


    if [[ $biv_only == true ]]
    then

    if [[ $midseptum == false ]]
    then
    cmd="/home/common/CARPentry_KCL_latest/bin/ekbatch "$mesh_path"/BiV_FEC_w5_h"$FEC_lvl"_retagged_"$suffix" "$init_files_path"/"$SA_folder"_BiV.RV_endo.apex"

    elif [[ $midseptum == true ]]
    then
    cmd="/home/common/CARPentry_KCL_latest/bin/ekbatch "$mesh_path"/BiV_FEC_w5_h"$FEC_lvl"_retagged_"$suffix" "$init_files_path"/"$SA_folder"_BiV.midseptum"

    fi
    echo $cmd
    eval $cmd

    elif [[ $biv_only == false ]]
    then
    cmd="/home/common/CARPentry_KCL_latest/bin/ekbatch "$mesh_path"/BiV_FEC_w5_h"$FEC_lvl"_retagged_"$suffix" "
    for PHI in AN AL LA IL IN
    do
    for Z in `seq 1 8`
    do
    cmd=$cmd$init_files_path"/"$SA_folder"_"$PHI"_"$Z","
    done;
    done;

    cmd_no_last_char=${cmd%?}

    echo $cmd_no_last_char
    eval $cmd_no_last_char

    fi



done
