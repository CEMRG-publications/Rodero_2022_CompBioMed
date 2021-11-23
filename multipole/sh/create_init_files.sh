#!/bin/bash

clear

which_cases="both"
biv_only=false
midseptum=false
all_cases="01 02 03 04 05 06 07 08 09 "`seq 10 20`
CV_l=0.5
CV_t=0.3
kFEC=5
SA_folder="default_HF"

function Usage(){
    echo "Script to create the init files for the multipole simulations."
    echo "Parameters:"
    echo "--SA_folder: Name/preffix of the simulations"
    echo "--CV_l: Conduction velocity in the longitudinal direction of the fibres."
    echo "--CV_t: Conduction velocity in the transverse direction of the fibres."
    echo "--kFEC: Anysotropy in the endocardium."
    echo "--which_cases: h/RR, HF or both for healthy/reverse remodelled, heart failure of both, respectively."
    echo "--biv_only: If true, only created for the BiV RV apex, false creates all the others. Default is true."
    echo "--midseptum: If true, it creates the init file for the midseptum electrode."
    echo "-h/--help: Parameters usage."
}


while [ "$1" != "" ]; do
    case $1 in
        --CV_l)                 shift
                                CV_l=$1
                                ;;
        --SA_folder)            shift
                                SA_folder=$1
                                ;;
        --CV_t)                 shift
                                CV_t=$1
                                ;;
        --kFEC)                 shift
                                kFEC=$1
                                ;;
        --which_cases)          shift
                                which_cases=$1
                                ;;
        --biv_only )            shift
                                biv_only=$1
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
cmd="/home/crg17/Desktop/scripts/multipole/sh/create_init_files.sh --which_cases RR --biv_only "$biv_only" --midseptum "$midseptum" --CV_l "$CV_l" --CV_t "$CV_t" --kFEC "$kFEC" --SA_folder "$SA_folder
echo $cmd
eval $cmd

cmd="/home/crg17/Desktop/scripts/multipole/sh/create_init_files.sh --which_cases HF --biv_only "$biv_only" --midseptum "$midseptum" --CV_l "$CV_l" --CV_t "$CV_t" --kFEC "$kFEC" --SA_folder "$SA_folder
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

common_path="/media/crg17/Seagate\ Backup\ Plus\ Drive/CT_cases"

for heart in $all_cases
    do
    init_files_path=$common_path"/"$which_cases"_case"$heart"/simulations/multipole/init_files"
    electrodes_path=$common_path"/"$which_cases"_case"$heart"/simulations/multipole/electrodes"

    cmd="mkdir -p "$init_files_path
    echo $cmd
    eval $cmd

    if [[ $biv_only == true ]]
    then

    if [[ $midseptum == false ]]
    then
    cmd="/home/crg17/Desktop/scripts/multipole/python/CARPtoEKBATCH.py
    --stimulus "$electrodes_path"/BiV.RV_endo.apex.vtx
    --output "$init_files_path"/"$SA_folder"_BiV.RV_endo.apex
    --CV_l "$CV_l" --CV_t "$CV_t" --kFEC "$kFEC
    echo $cmd
    eval $cmd

    elif [[ $midseptum == true ]]
    then
    cmd="/home/crg17/Desktop/scripts/multipole/python/CARPtoEKBATCH.py
    --stimulus "$electrodes_path"/BiV.midseptum.vtx
    --output "$init_files_path"/"$SA_folder"_BiV.midseptum
    --CV_l "$CV_l" --CV_t "$CV_t" --kFEC "$kFEC
    echo $cmd
    eval $cmd

    fi

    elif [[ $biv_only == false ]]
    then

    for PHI in AN AL LA IL IN
    do
    for Z in `seq 1 8`
    do
    cmd="/home/crg17/Desktop/scripts/multipole/python/CARPtoEKBATCH.py
    --stimulus "$electrodes_path"/"$PHI"_"$Z".vtx
    --output "$init_files_path"/"$SA_folder"_"$PHI"_"$Z"
    --CV_l "$CV_l" --CV_t "$CV_t" --kFEC "$kFEC

    echo $cmd
    eval $cmd

    done;
    done;

    fi



done
