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

cmd="/home/crg17/Desktop/scripts/multipole/sh/run_multipole_ekbatch.sh --SA_folder scar --biv_only false --suffix noPVTV_scar5mm --which_cases HF"
echo $cmd
eval $cmd


/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh "Simulations finished"