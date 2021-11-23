#!/bin/bash

# /home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh "Starting the creation of the init files"
# cmd="/home/crg17/Desktop/scripts/multipole/sh/create_init_files.sh --SA_folder scar --biv_only true"
# echo $cmd
# eval $cmd

# cmd="/home/crg17/Desktop/scripts/multipole/sh/create_init_files.sh --SA_folder scar --biv_only false"
# eval $cmd


# /home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh "init files created. Now select the appropriate mesh in the meshes folder!"



# /home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh "Running ekbatch simulations"
# cmd="/home/crg17/Desktop/scripts/multipole/sh/run_multipole_ekbatch.sh --SA_folder scar --biv_only true --suffix noPVTV_scar5mm --which_cases HF"
# echo $cmd
# eval $cmd

#cmd="/home/crg17/Desktop/scripts/multipole/sh/run_multipole_ekbatch.sh --SA_folder scar --biv_only false --suffix noPVTV_scar5mm --which_cases HF"
#echo $cmd
#eval $cmd
#
#
#/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh "Simulations finished"

###### Re do everything with healthy parameters
#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/create_init_files.sh --SA_folder RR_param --CV_l 0.8 --CV_t 0.23 --kFEC 7 --which_cases RR --biv_only true"
#
#echo $cmd
#eval $cmd
#
#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/create_init_files.sh --SA_folder RR_param --CV_l 0.8 --CV_t 0.23 --kFEC 7 --which_cases RR --biv_only false"
#
#echo $cmd
#eval $cmd
#
#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/create_init_files.sh --SA_folder RR_param --CV_l 0.8 --CV_t 0.23 --kFEC 7 --which_cases HF --biv_only true"
#
#echo $cmd
#eval $cmd
#
#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/create_init_files.sh --SA_folder RR_param --CV_l 0.8 --CV_t 0.23 --kFEC 7 --which_cases HF --biv_only false"
#
#echo $cmd
#eval $cmd
#
#/home/crg17/Desktop/KCL_projects/MPP/4chmodel/sh/sendmail.sh "Init files created"

#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/run_multipole_ekbatch.sh --SA_folder RR_param --FEC_lvl 70 --which_cases RR --biv_only false"
#echo $cmd
#eval $cmd
#
#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/run_multipole_ekbatch.sh --SA_folder RR_param --FEC_lvl 70 --which_cases RR --biv_only true"
#echo $cmd
#eval $cmd
#
#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/run_multipole_ekbatch.sh --SA_folder RR_param --FEC_lvl 70 --which_cases HF --biv_only false"
#echo $cmd
#eval $cmd
#
#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/run_multipole_ekbatch.sh --SA_folder RR_param --FEC_lvl 70 --which_cases HF --biv_only true"
#echo $cmd
#eval $cmd

##### Re-do everyhting with big_RR and small_HF

#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/create_init_files.sh --SA_folder big_RR --which_cases RR --biv_only true"
#
#echo $cmd
#eval $cmd
#
#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/create_init_files.sh --SA_folder big_RR --which_cases RR --biv_only false"
#
#echo $cmd
#eval $cmd
#
#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/create_init_files.sh --SA_folder small_HF --which_cases HF --biv_only true"
#
#echo $cmd
#eval $cmd
#
#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/create_init_files.sh --SA_folder small_HF --which_cases HF --biv_only false"
#
#echo $cmd
#eval $cmd


#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/run_multipole_ekbatch.sh --SA_folder big_RR --suffix big --which_cases RR --biv_only false"
#echo $cmd
#eval $cmd
#
#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/run_multipole_ekbatch.sh --SA_folder big_RR --suffix big --which_cases RR --biv_only true"
#echo $cmd
#eval $cmd

#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/run_multipole_ekbatch.sh --SA_folder small_HF --suffix small --which_cases HF --biv_only false"
#echo $cmd
#eval $cmd
#
#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/run_multipole_ekbatch.sh --SA_folder small_HF --suffix small --which_cases HF --biv_only true"
#echo $cmd
#eval $cmd

##### Re-do everyhting with scar of 6 mm thick

#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/create_init_files.sh --SA_folder scar_6mm --which_cases HF --biv_only true"
#
#echo $cmd
#eval $cmd
#
#cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/create_init_files.sh --SA_folder scar_6mm --which_cases HF --biv_only false"
#
#echo $cmd
#eval $cmd


cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/run_multipole_ekbatch.sh --SA_folder scar_6mm --suffix noPVTV_scar6mm --which_cases HF --biv_only false"
echo $cmd
eval $cmd

cmd="/home/crg17/Desktop/KCL_projects/MPP/multipole/sh/run_multipole_ekbatch.sh --SA_folder scar_6mm --suffix noPVTV_scar6mm --which_cases HF --biv_only true"
echo $cmd
eval $cmd